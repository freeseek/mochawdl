version development

## Copyright (c) 2021-2025 Giulio Genovese
##
## Version 2025-08-19
##
## Contact Giulio Genovese <giulio.genovese@gmail.com>
##
## This WDL workflow computes poligenic risk scores
##
## Cromwell version support
## - Successfully tested on v90
##
## Distributed under terms of the MIT License

struct Reference {
  File fasta_fai
  Int min_chr_len
}

workflow score {
  input {
    String sample_set_id
    String sample_header = "sample_id"
    String? region
    String? tag # GP, AP, HDS, DS, GT, AS
    String ext_string = "scores"
    File? colheaders_tsv_file
    String? summary_path
    Array[File] summary_files
    Array[File]? summary_idxs
    Array[Float]? q_score_thr
    File? covar_tsv_file

    String ref_name = "GRCh38"
    String? ref_path
    String? ref_fasta_fai
    Int? min_chr_len

    File impute_tsv_file # batch_id path chr1_imp_vcf chr1_imp_vcf_index chr2_imp_vcf chr2_imp_vcf_index ...
    String? impute_data_path
    File? samples_file
    String? exclude_str
    String? include_str
    String basic_bash_docker = "debian:stable-slim"
    String docker_repository = "us.gcr.io/mccarroll-mocha"
    String bcftools_docker = "bcftools:1.22-20250819"
    String r_mocha_docker = "r_mocha:1.22-20250819"
  }

  String docker_repository_with_sep = docker_repository + if docker_repository != "" && docker_repository == sub(docker_repository, "/$", "") then "/" else ""

  String summary_path_with_sep = select_first([summary_path, ""]) + if defined(summary_path) && select_first([summary_path]) == sub(select_first([summary_path]), "/$", "") then "/" else ""
  String ref_path_with_sep = select_first([ref_path, ""]) + if defined(ref_path) && select_first([ref_path]) == sub(select_first([ref_path]), "/$", "") then "/" else ""
  Reference ref = object {
    fasta_fai: ref_path_with_sep + select_first([ref_fasta_fai, if ref_name == "GRCh38" then "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai" else if ref_name == "GRCh37" then "human_g1k_v37.fasta.fai" else None]),
    min_chr_len: select_first([min_chr_len, 2000000]),
  }
  # call the relevant chromosome
  String? chr_string = if defined(region) then sub(sub(select_first([region]), ":.*$", ""), "^chr", "") else None

  # read table with batches information (scatter could be avoided if there was a tail() function)
  Array[Array[String]] impute_tsv = read_tsv(impute_tsv_file)
  Int n_batches = length(impute_tsv)-1
  scatter (idx in range(n_batches)) { Array[String] impute_tsv_rows = impute_tsv[(idx+1)] }
  Map[String, Array[String]] impute_tbl = as_map(zip(impute_tsv[0], transpose(impute_tsv_rows)))
  # check if path is in impute table (see http://github.com/openwdl/wdl/issues/305)
  Boolean is_path_in_impute_tbl = length(collect_by_key(zip(flatten([keys(impute_tbl),["path"]]),range(length(keys(impute_tbl))+1)))["path"])>1

  # compute data paths for each batch
  scatter (idx in range(n_batches)) {
    String impute_data_paths_with_sep = (if defined(impute_data_path) then sub(select_first([impute_data_path]), "/$", "") + "/" else "") +
                                        (if is_path_in_impute_tbl then sub(impute_tbl["path"][idx], "/$", "") + "/" else "")
  }

  Array[Array[String]] ref_fasta_fai_tbl = transpose(read_tsv(ref.fasta_fai))
  scatter (idx in range(length(ref_fasta_fai_tbl[0]))) {
    Int fai_len = ref_fasta_fai_tbl[1][idx]
    if (fai_len > ref.min_chr_len && ref_fasta_fai_tbl[0][idx] != "Y" && ref_fasta_fai_tbl[0][idx] != "chrY") {
      String chr_strings = sub(ref_fasta_fai_tbl[0][idx], "^chr", "")
    }
  }

  scatter (p in cross(range(n_batches), range(if defined(region) then 1 else length(select_all(chr_strings))))) {
    Int batch_idx = p.left
    String hdr = "chr" + (if defined(region) then select_first([chr_string]) else select_all(chr_strings)[p.right]) + "_imp_vcf"
    call vcf_score {
      input:
        vcf_file = impute_data_paths_with_sep[p.left] + impute_tbl[hdr][p.left],
        vcf_idx = impute_data_paths_with_sep[p.left] + impute_tbl[(hdr + "_index")][p.left],
        samples_file = samples_file,
        q_score_thr = q_score_thr,
        colheaders_tsv_file = colheaders_tsv_file,
        summary_files = prefix(summary_path_with_sep, summary_files),
        summary_idxs = if defined(summary_idxs) then prefix(summary_path_with_sep, select_first([summary_idxs])) else None,
        sample_header = sample_header,
        region = region,
        tag = tag,
        exclude_str = exclude_str,
        include_str = include_str,
        filebase = basename(basename(impute_data_paths_with_sep[p.left] + impute_tbl[hdr][p.left], ".bcf"), ".vcf.gz") + "." + ext_string,
        docker = docker_repository_with_sep + bcftools_docker
    }
  }

  call score_summary {
    input:
      score_files = vcf_score.file,
      covar_tsv_file = covar_tsv_file,
      sample_header = sample_header,
      filebase = sample_set_id,
      ext_string = ext_string,
      docker = if defined(covar_tsv_file) then docker_repository_with_sep + r_mocha_docker else basic_bash_docker
  }

  output {
    File score_tsv_file = score_summary.file
    File? adj_score_tsv_file = score_summary.adj_file
  }

  meta {
    author: "Giulio Genovese"
    email: "giulio.genovese@gmail.com"
    description: "See the [MoChA](http://github.com/freeseek/mocha) website for more information"
  }
}

task vcf_score {
  input {
    File vcf_file
    File vcf_idx
    File? samples_file
    File? colheaders_tsv_file
    Array[File]+ summary_files
    Array[File]? summary_idxs
    Array[Float]? q_score_thr
    String? tag
    String? sample_header
    String? region
    String? exclude_str
    String? include_str
    String filebase

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Float summary_size = size(summary_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + vcf_size + summary_size)])
  Float memory = select_first([memory_override, 3.5 + if defined(summary_idxs) then 0 else summary_size])

  command <<<
    set -euo pipefail
    summary_files=~{write_lines(summary_files)}
    summary_idxs=~{if defined(summary_idxs) then write_lines(select_first([summary_idxs])) else ""}
    echo "~{sep("\n", select_all([vcf_file, vcf_idx, samples_file, colheaders_tsv_file]))}" | \
      cat - $summary_files ~{if defined(summary_idxs) then "$summary_idxs" else ""}| tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $summary_files
    bcftools +score \
      --summaries $summary_files \
      ~{if defined(tag) then "--use \"" + tag + "\"" else ""} \
      ~{if defined(q_score_thr) then "--q-score-thr " + sep(",", select_first([q_score_thr])) else ""} \
      --output "~{filebase}.tsv" \
      ~{if defined(sample_header) then "--sample-header \"" + sample_header + "\"" else ""} \
      ~{if defined(exclude_str) then "--exclude '" + exclude_str + "'" else ""} \
      ~{if defined(include_str) then "--include '" + include_str + "'" else ""} \
      ~{if defined(region) then "--regions \"" + select_first([region]) + "\"" else ""} \
      ~{if defined(samples_file) then "--samples-file \"" + basename(select_first([samples_file])) + "\" \\\n" +
      "  --force-samples" else ""} \
      ~{if defined(colheaders_tsv_file) then "--columns-file \"" + basename(select_first([colheaders_tsv_file])) + "\"" else ""} \
      "~{basename(vcf_file)}"
    echo "~{sep("\n", select_all([vcf_file, vcf_idx, samples_file, colheaders_tsv_file]))}" | \
      cat - $summary_files ~{if defined(summary_idxs) then "$summary_idxs" else ""}| sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File file = filebase + ".tsv"
  }

  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    docker: docker
    preemptible: preemptible
    maxRetries: maxRetries
  }
}

task score_summary {
  input {
    Array[File]+ score_files
    File? covar_tsv_file
    String sample_header
    String filebase
    String ext_string

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float score_size = size(score_files[0], "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + (length(score_files) + 1) * score_size)])
  Float memory = select_first([memory_override, 3.5 + 4.0 * score_size])

  command <<<
    set -euo pipefail
    score_files=~{write_lines(score_files)}
    cat - $score_files | tr '\n' '\0' | xargs -0 mv -t .
    ~{if defined(covar_tsv_file) then "mv \"" + select_first([covar_tsv_file]) + "\" ." else ""}
    cat $score_files | sed 's/^.*\///' | tr '\n' '\0' | xargs -0 cat | \
      awk 'NR==1 {sample_header=$1}
      {if ($1==sample_header) {
        for (i=2; i<=NF; i++) {
          if (!($i in score_id)) {
            col_count++;
            cols[col_count]=$i;
            score_id[$i]++;
          }
          f[i] = $i;
        }
      } else {
        if (!($1 in sample_id)) {
          row_count++;
          rows[row_count]=$1;
          sample_id[$1]++;
        }
        for (i=2; i<=NF; i++)
          v[f[i]"~"$1]+=$i;
      }}
      END {printf "%s",sample_header;
        for (i=1; i<=col_count; i++)
          printf "\t%s",cols[i];
        printf "\n";
        for(j=1; j<=row_count; j++) {
          printf "%s",rows[j];
          for (i=1; i<=col_count; i++)
            printf "\t%f",v[cols[i]"~"rows[j]];
          printf "\n";
        }
      }' > "~{filebase}.~{ext_string}.tsv"
    ~{if defined(covar_tsv_file) then
    "R --vanilla <<CODE\n" +
    "library(data.table)\n" +
    "df_scores <- fread('" + filebase + "." + ext_string + ".tsv', sep = \"\\t\", header = TRUE, data.table = FALSE)\n" +
    "df_covars <- fread('" + basename(select_first([covar_tsv_file])) + "', sep = \"\\t\", header = TRUE, data.table = FALSE)\n" +
    "df <- merge(df_scores, df_covars, by = '" + sample_header + "')\n" +
    "df_adj <- data.frame(" + sample_header + " = df[, '" + sample_header + "'])\n" +
    "scores <- names(df_scores)\n" +
    "scores <- scores[scores != '" + sample_header + "']\n" +
    "covars <- names(df_covars)\n" +
    "covars <- covars[covars != '" + sample_header + "']\n" +
    "bt <- rawToChar(as.raw(96))\n" +
    "for (score in scores) {\n" +
    "  formula <- paste0(bt, score, bt, ' ~ ', bt, paste(covars, collapse = paste(bt, '+', bt)), bt)\n" +
    "  fit <- lm(formula, df)\n" +
    "  df_adj[, paste0('adj_', score)] <- df[, score] - predict(fit, df)\n" +
    "  df_adj[, paste0('adj_', score)] <- df_adj[, paste0('adj_', score)] / sd(df_adj[, paste0('adj_', score)])\n" +
    "}\n" +
    "write.table(df_adj, '" + filebase + ".adj_" + ext_string + ".tsv', sep = '\t', quote = FALSE, row.names = FALSE)\n" +
    "CODE" else ""}
    cat - $score_files | sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
    ~{if defined(covar_tsv_file) then "rm \"" + basename(select_first([covar_tsv_file])) + "\"" else ""}
  >>>

  output {
    File file = filebase + "." + ext_string + ".tsv"
    File? adj_file = if defined(covar_tsv_file) then filebase + ".adj_" + ext_string + ".tsv" else None
  }

  runtime {
    docker: docker
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GiB"
    preemptible: preemptible
    maxRetries: maxRetries
  }
}
