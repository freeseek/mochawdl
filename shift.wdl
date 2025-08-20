version development

## Copyright (c) 2021-2025 Giulio Genovese
##
## Version 2025-08-19
##
## Contact Giulio Genovese <giulio.genovese@gmail.com>
##
## This WDL workflow runs allelic shift imbalance analysis in a given region
##
## Cromwell version support
## - Successfully tested on v90
##
## Distributed under terms of the MIT License

struct Reference {
  File fasta
  File fasta_fai
  String n_x_chr
  File? cyto_file
  String? chr_prefix
  Array[Int] len
  Array[Int] pcen
  Array[Int] qcen
  File? gff3_file # http://ftp.ensembl.org/pub/current_gff3/homo_sapiens/
  File? rsid_vcf_file # http://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/
  File? rsid_vcf_idx
}

workflow shift {
  input {
    String sample_set_id
    File? keep_samples_file
    File? remove_samples_file
    File pheno_tsv_file
    String as_id = "AS"
    String ext_string = "as"
    String? pop

    String ref_name = "GRCh38"
    String? ref_fasta
    String? ref_fasta_fai
    String? ref_path
    String? chr_prefix
    Int? n_x_chr
    String? cyto_file
    String? gff3_file
    String? rsid_vcf_file
    String? rsid_vcf_idx

    File impute_tsv_file # batch_id path chr1_imp_vcf chr1_imp_vcf_index chr2_imp_vcf chr2_imp_vcf_index ...
    String? impute_data_path
    Boolean fisher_exact = true
    Boolean drop_genotypes = true
    Boolean phred_score = true
    Boolean plot = true
    String basic_bash_docker = "debian:stable-slim"
    String docker_repository = "us.gcr.io/mccarroll-mocha"
    String bcftools_docker = "bcftools:1.22-20250819"
    String r_mocha_docker = "r_mocha:1.22-20250819"
  }

  String docker_repository_with_sep = docker_repository + if docker_repository != "" && docker_repository == sub(docker_repository, "/$", "") then "/" else ""

  String ref_path_with_sep = select_first([ref_path, ""]) + if defined(ref_path) && select_first([ref_path]) == sub(select_first([ref_path]), "/$", "") then "/" else ""
  Reference ref = object {
    fasta: if defined(ref_fasta) then ref_path_with_sep + select_first([ref_fasta]) else if ref_name == "GRCh38" then ref_path_with_sep + "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" else if ref_name == "GRCh37" then ref_path_with_sep + "human_g1k_v37.fasta" else None,
    fasta_fai: if defined(ref_fasta_fai) then ref_path_with_sep + select_first([ref_fasta_fai]) else if ref_name == "GRCh38" then ref_path_with_sep + "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai" else if ref_name == "GRCh37" then ref_path_with_sep + "human_g1k_v37.fasta.fai" else None,
    n_x_chr: select_first([n_x_chr, 23]),
    cyto_file: if defined(ref_path) || defined(cyto_file) then ref_path_with_sep + select_first([cyto_file, "cytoBand.txt.gz"]) else None,
    chr_prefix: if ref_name == "GRCh38" then "chr" else if ref_name == "GRCh37" then "" else chr_prefix,
    len: if ref_name == "GRCh38" then
      [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895]
    else if ref_name == "GRCh37" then
      [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560]
    else None,
    pcen: if ref_name == "GRCh38" then
      [122026459, 92188145, 90772458, 49712061, 46485900, 58553888, 58169653, 44033744, 43389635, 39686682, 51078348, 34769407, 16000000, 16000000, 17083673, 36311158, 22813679, 15460899, 24498980, 26436232, 10864560, 12954788, 58605579]
    else if ref_name == "GRCh37" then
      [121535434, 92326171, 90504854, 49660117, 46405641, 58830166, 58054331, 43838887, 47367679, 39254935, 51644205, 34856694, 16000000, 16000000, 17000000, 35335801, 22263006, 15460898, 24681782, 26369569, 11288129, 13000000, 58632012]
    else None,
    qcen: if ref_name == "GRCh38" then
      [124932724, 94090557, 93655574, 51743951, 50059807, 59829934, 61528020, 45877265, 45518558, 41593521, 54425074, 37185252, 18051248, 18173523, 19725254, 38265669, 26616164, 20861206, 27190874, 30038348, 12915808, 15054318, 62412542]
    else if ref_name == "GRCh37" then
      [124535434, 95326171, 93504854, 52660117, 49405641, 61830166, 61054331, 46838887, 50367679, 42254935, 54644205, 37856694, 19000000, 19000000, 20000000, 38335801, 25263006, 18460898, 27681782, 29369569, 14288129, 16000000, 61632012]
    else None,
    gff3_file: if defined(gff3_file) then ref_path_with_sep + select_first([gff3_file]) else None,
    rsid_vcf_file: if defined(rsid_vcf_file) then ref_path_with_sep + select_first([rsid_vcf_file]) else None,
    rsid_vcf_idx: if defined(rsid_vcf_idx) then ref_path_with_sep + select_first([rsid_vcf_idx]) else None
  }

  # read table with batches information (scatter could be avoided if there was a tail() function)
  Array[Array[String]] impute_tsv = read_tsv(impute_tsv_file)
  Int n_batches = length(impute_tsv)-1
  scatter (idx in range(n_batches)) { Array[String] impute_tsv_rows = impute_tsv[(idx+1)] }
  Map[String, Array[String]] impute_tbl = as_map(zip(impute_tsv[0], transpose(impute_tsv_rows)))
  # check if path is in impute table (see http://github.com/openwdl/wdl/issues/305)
  Boolean is_path_in_impute_tbl = length(collect_by_key(zip(flatten([keys(impute_tbl),["path"]]),range(length(keys(impute_tbl))+1)))["path"])>1

  # compute data paths for each batch
  scatter (idx in range(n_batches)) {
    String impute_data_paths = select_first([impute_data_path, if is_path_in_impute_tbl then impute_tbl["path"][idx] else ""])
    String impute_data_paths_with_sep = impute_data_paths + (if impute_data_paths == "" || sub(impute_data_paths, "/$", "") != impute_data_paths then "" else "/")
  }

  call lst_header { input: pheno_tsv_file = pheno_tsv_file, docker = basic_bash_docker }
  # compute phenotype regions to test
  scatter (idx in range(length(lst_header.phenos))) {
    String region_name = sub(lst_header.phenos[idx], "_.*$", "")
    String chr_string = sub(sub(region_name, "[pq]*$", ""), "Y", "X")
    if (sub(chr_string, "[0-9X]+", "") == "") {
      Int pheno_idx = idx
      String pheno_chr = chr_string
      String pheno_name = lst_header.phenos[idx]
      String vcf_file_suffix = pheno_name + "." + ext_string + ".bcf" # http://github.com/broadinstitute/cromwell/issues/5549
      String vcf_idx_suffix = pheno_name + "." + ext_string + ".bcf.csi" # http://github.com/broadinstitute/cromwell/issues/5549
      Int chr_idx = sub(chr_string, "^X$", ref.n_x_chr)
      String arm = sub(region_name, "^[0-9XY]*", "")
      String pheno_regions = ref.chr_prefix + chr_string + if arm == "p" then ":1-" + ref.pcen[(chr_idx - 1)] else if arm == "q" then ":" + ref.qcen[(chr_idx - 1)] + "-" + ref.len[(chr_idx - 1)] else ""
    }
  }

  # generate list of all chromosomes to be processed
  Map[String, Array[String]] chr2pheno_idx = collect_by_key(zip(select_all(pheno_chr), select_all(pheno_idx)))
  Map[String, Array[String]] chr2pheno_names = collect_by_key(zip(select_all(pheno_chr), select_all(pheno_name)))
  Map[String, Array[String]] chr2regions = collect_by_key(zip(select_all(pheno_chr), select_all(pheno_regions)))
  Map[String, Array[String]] chr2vcf_file_suffix = collect_by_key(zip(select_all(pheno_chr), select_all(vcf_file_suffix))) # http://github.com/broadinstitute/cromwell/issues/5549
  Map[String, Array[String]] chr2vcf_idx_suffix = collect_by_key(zip(select_all(pheno_chr), select_all(vcf_idx_suffix))) # http://github.com/broadinstitute/cromwell/issues/5549
  scatter (p in cross(range(n_batches), keys(chr2pheno_names))) {
    Array[Int] cross_idx = chr2pheno_idx[p.right]
    call vcf_summary {
      input:
        vcf_file = impute_data_paths_with_sep[p.left] + impute_tbl[("chr" + p.right + "_imp_vcf")][p.left],
        vcf_idx = impute_data_paths_with_sep[p.left] + impute_tbl[("chr" + p.right + "_imp_vcf_index")][p.left],
        pheno_names = chr2pheno_names[p.right],
        regions = chr2regions[p.right],
        vcf_file_suffix = chr2vcf_file_suffix[p.right], # http://github.com/broadinstitute/cromwell/issues/5549
        vcf_idx_suffix = chr2vcf_idx_suffix[p.right], # http://github.com/broadinstitute/cromwell/issues/5549
        keep_samples_file = keep_samples_file,
        remove_samples_file = remove_samples_file,
        pheno_tsv_file = pheno_tsv_file,
        fisher_exact = fisher_exact,
        filebase = sample_set_id + "." + p.left,
        as_id = as_id,
        ext_string = ext_string,
        drop_genotypes = drop_genotypes,
        docker = docker_repository_with_sep + bcftools_docker
    }
  }

  scatter (idx in select_all(pheno_idx)) {
    Map[Int, Array[File]] idx2vcf_files = collect_by_key(zip(flatten(cross_idx), flatten(vcf_summary.as_vcf_files)))
    Map[Int, Array[File]] idx2vcf_idxs = collect_by_key(zip(flatten(cross_idx), flatten(vcf_summary.as_vcf_idxs)))
    call vcf_merge {
      input:
        vcf_files = idx2vcf_files[idx],
        vcf_idxs = idx2vcf_idxs[idx],
        ref_fasta = if defined(ref.gff3_file) then ref.fasta else None,
        fasta_fai = ref.fasta_fai,
        gff3_file = ref.gff3_file,
        rsid_vcf_file = ref.rsid_vcf_file,
        rsid_vcf_idx = ref.rsid_vcf_idx,
        region = select_first([pheno_regions[idx]]),
        fisher_exact = fisher_exact,
        as_id = as_id,
        filebase = sample_set_id + "." + ext_string + "_" + lst_header.phenos[idx] + if defined(pop) then "." + select_first([pop]) else "",
        pheno_name = ext_string + "_" + lst_header.phenos[idx] + if defined(pop) then "." + select_first([pop]) else "",
        phred_score = phred_score,
        docker = docker_repository_with_sep + bcftools_docker
    }

    if (plot) {
      if (vcf_merge.n > 0) {
        call assoc_plot {
          input:
            vcf_file = vcf_merge.gwas_vcf_file,
            vcf_idx = vcf_merge.gwas_vcf_idx,
            region = select_first([pheno_regions[idx]]),
            cyto_file = ref.cyto_file,
            csq = defined(ref.gff3_file),
            filebase = sample_set_id + "." + ext_string + "_" + lst_header.phenos[idx] + if defined(pop) then "." + select_first([pop]) else "",
            docker = docker_repository_with_sep + r_mocha_docker
        }
      }
    }
  }

  output {
    Array[File] gwas_vcf_file = vcf_merge.gwas_vcf_file
    Array[File] gwas_vcf_idx = vcf_merge.gwas_vcf_idx
    Array[File]? png_files = if plot then select_all(assoc_plot.png_file) else None
  }

  meta {
    author: "Giulio Genovese"
    email: "giulio.genovese@gmail.com"
    description: "See the [MoChA](http://github.com/freeseek/mocha) website for more information"
  }
}

task lst_header {
  input {
    File pheno_tsv_file
    String space_character = '_'

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mv "~{pheno_tsv_file}" .
    head -n1 "~{basename(pheno_tsv_file)}" | cut -f2- | tr '\t ' '\n~{space_character}'
    rm "~{basename(pheno_tsv_file)}"
  >>>

  output {
    Array[String] phenos = read_lines(stdout())
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

# the command requires BCFtools 1.14 due to bug http://github.com/samtools/bcftools/issues/1566
task vcf_summary {
  input {
    File vcf_file
    File vcf_idx
    Array[String]+ pheno_names
    Array[String]+ regions
    Array[String]+ vcf_file_suffix # suffix array passed due to bug http://github.com/broadinstitute/cromwell/issues/5549
    Array[String]+ vcf_idx_suffix # suffix array passed due to bug http://github.com/broadinstitute/cromwell/issues/5549
    File? keep_samples_file
    File? remove_samples_file
    File pheno_tsv_file
    Boolean fisher_exact
    String filebase
    String as_id
    String ext_string
    Boolean drop_genotypes = true

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size)])

  command <<<
    set -euo pipefail
    echo "~{sep("\n", select_all([vcf_file, vcf_idx, keep_samples_file, remove_samples_file, pheno_tsv_file]))}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    pheno_names=~{write_lines(pheno_names)}
    regions=~{write_lines(regions)}
    paste $pheno_names $regions | while IFS=$'\t' read pheno_name region; do
      awk -F"\t" -v pheno_name=$pheno_name 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
        NR>1 && ($(f[pheno_name])==0 || $(f[pheno_name])==1) {print $(f["sample_id"])"\t"$(f[pheno_name])}' \
        "~{basename(pheno_tsv_file)}" > "~{filebase}.pheno.tsv"
      cut -f1 "~{filebase}.pheno.tsv"~{if defined(keep_samples_file) then " | \\\n" +
      "awk -F\"\\t\" 'NR==FNR {x[$1]++} NR>FNR && $1 in x' \"" + basename(select_first([keep_samples_file])) + "\" -"
      else ""}~{if defined(remove_samples_file) then " | \\\n" +
      "awk -F\"\\t\" 'NR==FNR {x[$1]++} NR>FNR && !($1 in x)' \"" + basename(select_first([remove_samples_file])) + "\" -"
      else ""} > "~{filebase}.samples.lines"
      ~{if fisher_exact then "awk -F\"\\t\" '$2==0 {print $1}'  \"" + filebase + ".pheno.tsv\" > \"" + filebase + ".controls.lines\"" else ""}
      awk -F"\t" '$2==1 {print $1}' "~{filebase}.pheno.tsv" > "~{filebase}.cases.lines"
      bcftools +mochatools \
        --no-version \
        --output-type u \
        --regions $region \
        "~{basename(vcf_file)}" \
        -- --tags MACH \
        --samples-file "~{filebase}.samples.lines" \
        --force-samples | \
      bcftools annotate \
        --no-version \
        --output-type u \
        --remove ID,QUAL,FILTER,^INFO/MACH,^FMT/GT,FMT/~{as_id}~{if fisher_exact then " | \\\n" +
      "bcftools +contrast \\\n" +
      "  --output-type u \\\n" +
      "  --annots NASSOC \\\n" +
      "  --control-samples \"" + filebase + ".controls.lines\" \\\n" +
      "  --case-samples \"" + filebase + ".cases.lines\" \\\n" +
      "  --force-samples"
      else ""} | \
      bcftools +mochatools \
        --no-version \
        --output-type u \
        -- --summary ~{as_id} \
        --samples-file "~{filebase}.cases.lines" \
        --force-samples \
        ~{if drop_genotypes then "--drop-genotypes" else ""} | \
      bcftools view --no-version --output "~{filebase}.$pheno_name.~{ext_string}.bcf" --output-type b --include 'sum(~{as_id})>0' --write-index
      rm "~{filebase}.pheno.tsv" "~{filebase}.samples.lines" ~{if fisher_exact then "\"" + filebase + ".controls.lines\" " else ""}"~{filebase}.cases.lines"
    done
    echo "~{sep("\n", select_all([vcf_file, vcf_idx, keep_samples_file, remove_samples_file, pheno_tsv_file]))}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    Array[File] as_vcf_files = prefix(filebase + ".", vcf_file_suffix)
    Array[File] as_vcf_idxs = prefix(filebase + ".", vcf_idx_suffix)
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

# bcftools annotate needs --regions due to bug http://github.com/samtools/bcftools/issues/1199
# see http://github.com/MRCIEU/gwas-vcf-specification for VCF output
# see Marchini, J., Howie, B. Genotype imputation for genome-wide association studies. Nat Rev Genet 11, 499–511 (2010). http://doi.org/10.1038/nrg2796
task vcf_merge {
  input {
    Array[File]+ vcf_files
    Array[File]+ vcf_idxs
    File? ref_fasta
    File fasta_fai
    File? gff3_file
    File? rsid_vcf_file
    File? rsid_vcf_idx
    String region
    Boolean fisher_exact
    String as_id
    String filebase
    String pheno_name
    Boolean phred_score = true

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_files, "GiB")
  Float ref_size = size(ref_fasta, "GiB")
  Float dbsnp_size = size(rsid_vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size + ref_size + dbsnp_size)])

  command <<<
    set -euo pipefail
    vcf_files=~{write_lines(vcf_files)}
    vcf_idxs=~{write_lines(vcf_idxs)}
    echo "~{sep("\n", select_all([ref_fasta, fasta_fai, gff3_file, rsid_vcf_file, rsid_vcf_idx]))}" | \
      cat - $vcf_files $vcf_idxs | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $vcf_files $vcf_idxs
    bcftools merge \
      --no-version \
      --output-type u \
      --force-single \
      --info-rules MACH:sum,~{if fisher_exact then "NASSOC:sum," else ""}~{as_id}:sum \
      --file-list $vcf_files \
      --merge none | \
    bcftools reheader \
      --fai "~{basename(fasta_fai)}" \
      --temp-prefix ./bcftools. | \
    ~{if fisher_exact then
      "bcftools +mochatools \\\n" +
      "  --no-version \\\n" +
      "  --output-type u \\\n" +
      "  -- --test NASSOC \\\n" +
      (if phred_score then "  --phred" else "") + " |"
    else ""} \
    bcftools +mochatools \
      --no-version \
      --output-type ~{if defined(gff3_file) then "u" else "b"} \
      -- --test ~{as_id} \
      ~{if phred_score then "--phred" else ""} \
    ~{if defined(gff3_file) then "  | \\\n" +
    "bcftools csq \\\n" +
    "  --no-version \\\n" +
    "  --output-type b \\\n" +
    "  --fasta-ref \"" + basename(select_first([ref_fasta])) + "\" \\\n" +
    "  --gff-annot \"" + basename(select_first([gff3_file])) + "\" \\\n" +
    "  --trim-protein-seq 1 \\\n" +
    "  --custom-tag CSQ \\\n" +
    "  --local-csq \\\n" +
    "  --ncsq 64 \\\n" +
    "  --samples -" else ""} \
      --output "~{filebase + if defined(rsid_vcf_file) then ".tmp" else ""}.bcf" \
      --write-index
    ~{if defined(rsid_vcf_file) then
      "bcftools annotate \\\n" +
      "  --no-version \\\n" +
      "  --annotations \"" + basename(select_first([rsid_vcf_file])) + "\" \\\n" +
      "  --columns ID \\\n" +
      "  --output \"" + filebase + ".bcf\" \\\n" +
      "  --output-type b \\\n" +
      "  --regions \"" + region + "\" \\\n" +
      "  \"" + filebase + ".tmp.bcf\" \\\n" +
      "  --write-index\n" +
      "  rm \"" + filebase + ".tmp.bcf\" \"" + filebase + ".tmp.bcf.csi\""
      else ""}
    (echo "CHROM GENPOS ALLELE0 ALLELE1 N INFO N_CASES BETA SE LOG10P A1FREQ AC_ALLELE1";
    bcftools query --format "%CHROM\t%POS\t%REF\t%ALT\t%MACH{0}\t%MACH{1}\t%MACH{2}\t%AS{0}\t%AS{1}\t%pbinom_AS\n" "~{filebase}.bcf" | \
    awk -F"\t" '{if ($6*($5-$6/2)==0) si=1; else si=($5*$7/$6-$6)/($5-$6/2);
      nc=$8+$9; es=log(($9+.5)/($8+.5)); se=sqrt(1/($8+.5)+1/($9+.5)); lp=$10/10; af=$6/$5/2; ac=$6;
      print $1,$2,$3,$4,$5,si,nc,es,se,lp,af,ac}') | \
    bcftools +munge --no-version -c REGENIE --fai "~{basename(fasta_fai)}" -s ~{pheno_name} | \
    bcftools merge --no-version --merge none --no-index --output "~{filebase}.gwas.bcf" --output-type b "~{filebase}.bcf" --write-index -
    bcftools query -f "\n" -i 'sum(~{as_id})>1' "~{filebase}.gwas.bcf" | wc -l
    rm "~{filebase}.bcf" "~{filebase}.bcf.csi"
    echo "~{sep("\n", select_all([ref_fasta, fasta_fai, gff3_file, rsid_vcf_file, rsid_vcf_idx]))}" | \
      cat - $vcf_files $vcf_idxs | sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    Int n = read_int(stdout())
    File gwas_vcf_file = filebase + ".gwas.bcf"
    File gwas_vcf_idx = filebase + ".gwas.bcf.csi"
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

task assoc_plot {
  input {
    File vcf_file
    File vcf_idx
    String region
    File? cyto_file
    Boolean csq = false
    String filebase

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + vcf_size)])

  command <<<
    set -euo pipefail
    mv "~{vcf_file}" .
    mv "~{vcf_idx}" .
    ~{if defined(cyto_file) then "mv \"" + select_first([cyto_file]) + "\" ." else ""}
    assoc_plot.R \
      ~{if defined(cyto_file) then "--cytoband \"" + basename(select_first([cyto_file])) + "\"" else ""} \
      --vcf "~{basename(vcf_file)}" \
      --as \
      ~{if csq then "--csq" else ""} \
      --region ~{region} \
      --png "~{filebase}.png"
    rm "~{basename(vcf_file)}"
    rm "~{basename(vcf_idx)}"
    ~{if defined(cyto_file) then "rm \"" + basename(select_first([cyto_file])) + "\"" else ""}
  >>>

  output {
    File png_file = filebase + '.png'
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
