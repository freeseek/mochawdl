version development

## Copyright (c) 2021-2024 Giulio Genovese
##
## Version 2024-05-05
##
## Contact Giulio Genovese <giulio.genovese@gmail.com>
##
## This WDL workflow runs association analyses with REGENIE and PLINK2
##
## Cromwell version support
## - Successfully tested on v86
##
## Distributed under terms of the MIT License

struct Reference {
  File? fasta
  File fasta_fai
  Int min_chr_len
  Int n_x_chr
  Int? par_bp1
  Int? par_bp2
  File? cyto_file
  File genetic_map_file
  String? pca_exclusion_regions
  File? gff3_file # https://ftp.ensembl.org/pub/current_gff3/homo_sapiens/
  File? rsid_vcf_file # https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/
  File? rsid_vcf_idx
}

workflow assoc {
  input {
    String sample_set_id
    String? sex_specific # male female
    Float max_win_size_cm_step2 = 20.0
    File sample_tsv_file
    File? keep_samples_file
    File? remove_samples_file
    Int min_mac_step1 = 10
    Float min_maf_step1 = 0.01
    Int min_mac_step2 = 10
    Float? min_info_step2
    File? covar_tsv_file
    File? pheno_tsv_file
    String? pop
    String dosage_field = "DS"
    String space_character = "_"
    Boolean binary = true
    Int min_case_count = 20
    Int min_sex_count = 20
    Int bsize_step1 = 1000
    Int bsize_step2 = 400
    Int max_vif = 2000
    Float max_corr = 0.9999
    Float cis_plot_min_af = 0.01
    Boolean loocv = true
    String? regenie_step0_extra_args
    String? regenie_step1_extra_args
    String? regenie_step2_extra_args
    String? plink_extra_args
    Boolean step1 = true
    Boolean pgt_output = false
    Boolean pca = false
    Boolean step2 = true
    Boolean cis = false
    Boolean plot = true
    Int pca_ndim = 20
    Int pca_cpus = 2
    File? input_loco_lst
    String? input_loco_path
    File? input_firth_lst
    String? input_firth_path

    String ref_name = "GRCh38"
    String? ref_path
    String? ref_fasta
    String? ref_fasta_fai
    Int? min_chr_len
    Int? n_x_chr
    Int? par_bp1
    Int? par_bp2
    String? cyto_file
    String? genetic_map_file
    String? pca_exclusion_regions
    String? gff3_file
    String? rsid_vcf_file
    String? rsid_vcf_idx

    File? mocha_tsv_file # batch_id path pgt_vcf pgt_vcf_index
    String? mocha_data_path
    File? impute_tsv_file # batch_id path chr1_imp_vcf chr1_imp_vcf_index chr2_imp_vcf chr2_imp_vcf_index ...
    String? impute_data_path
    String basic_bash_docker = "debian:stable-slim"
    String pandas_docker = "amancevice/pandas:slim"
    String docker_repository = "us.gcr.io/mccarroll-mocha"
    String bcftools_docker = "bcftools:1.20-20240505"
    String regenie_docker = "regenie:1.20-20240505"
    String r_mocha_docker = "r_mocha:1.20-20240505"
  }

  String docker_repository_with_sep = docker_repository + if docker_repository != "" && docker_repository == sub(docker_repository, "/$", "") then "/" else ""

  String filebase = sub(sample_set_id, "[ \t]", "_") + if defined(sex_specific) then space_character + select_first([sex_specific]) else ""

  String ref_path_with_sep = select_first([ref_path, ""]) + if defined(ref_path) && select_first([ref_path]) == sub(select_first([ref_path]), "/$", "") then "/" else ""
  Reference ref = object {
    fasta: if defined(ref_fasta) then ref_path_with_sep + select_first([ref_fasta]) else if ref_name == "GRCh38" then ref_path_with_sep + "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" else if ref_name == "GRCh37" then ref_path_with_sep + "human_g1k_v37.fasta" else None,
    fasta_fai: if defined(ref_fasta_fai) then ref_path_with_sep + select_first([ref_fasta_fai]) else if ref_name == "GRCh38" then ref_path_with_sep + "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai" else if ref_name == "GRCh37" then ref_path_with_sep + "human_g1k_v37.fasta.fai" else None,
    min_chr_len: select_first([min_chr_len, 3000000]),
    n_x_chr: select_first([n_x_chr, 23]),
    par_bp1: if defined(par_bp1) then select_first([par_bp1]) else if ref_name == "GRCh38" then 2781479 else if ref_name == "GRCh37" then 2699520 else None,
    par_bp2: if defined(par_bp2) then select_first([par_bp2]) else if ref_name == "GRCh38" then 155701383 else if ref_name == "GRCh37" then 154931044 else None,
    cyto_file: if defined(ref_path) || defined(cyto_file) then ref_path_with_sep + select_first([cyto_file, "cytoBand.txt.gz"]) else None,
    genetic_map_file: if defined(genetic_map_file) then ref_path_with_sep + select_first([genetic_map_file]) else if ref_name == "GRCh38" then ref_path_with_sep + "genetic_map_hg38_withX.txt.gz" else if ref_name == "GRCh37" then ref_path_with_sep + "genetic_map_hg19_withX.txt.gz" else None,
    pca_exclusion_regions: if defined(pca_exclusion_regions) then pca_exclusion_regions else if ref_name == "GRCh38" then "5:43999898-52204166,6:24999772-33532223,8:8142478-12142491,11:44978449-57232526" else if ref_name == "GRCh37" then "5:44000000-51500000,6:25000000-33500000,8:8000000-12000000,11:45000000-57000000" else None,
    gff3_file: if defined(gff3_file) then ref_path_with_sep + select_first([gff3_file]) else None,
    rsid_vcf_file: if defined(rsid_vcf_file) then ref_path_with_sep + select_first([rsid_vcf_file]) else None,
    rsid_vcf_idx: if defined(rsid_vcf_idx) then ref_path_with_sep + select_first([rsid_vcf_idx]) else None
  }

  Array[Array[String]] ref_fasta_fai_tbl = transpose(read_tsv(ref.fasta_fai))
  scatter (idx in range(length(ref_fasta_fai_tbl[0]))) {
    Int fai_len = ref_fasta_fai_tbl[1][idx]
    if (fai_len > ref.min_chr_len && ref_fasta_fai_tbl[0][idx] != "Y" && ref_fasta_fai_tbl[0][idx] != "chrY") {
      String chrs = ref_fasta_fai_tbl[0][idx]
      Int lens = ref_fasta_fai_tbl[1][idx]
    }
  }

  if (defined(pheno_tsv_file)) {
    call prune_file {
      input:
        sex_specific = sex_specific,
        sample_tsv_file = sample_tsv_file,
        keep_samples_file = keep_samples_file,
        remove_samples_file = remove_samples_file,
        covar_tsv_file = covar_tsv_file,
        pheno_tsv_file = select_first([pheno_tsv_file]),
        space_character = space_character,
        binary = binary,
        min_case_count = min_case_count,
        min_sex_count = min_sex_count,
        filebase = filebase,
        docker = basic_bash_docker
    }
  }

  # REGENIE step 1

  if (step1 || pca) {
    call ref_scatter as ref_scatter_step1 {
      input:
        chrs = select_all(chrs),
        lens = select_all(lens),
        genetic_map_file = ref.genetic_map_file,
        max_win_size_cm = 300.0, # until regenie updates, step 1 cannot be parallelized beyond the 23 chromosomes
        overlap_size_cm = 0.0,
        genetic_map_order = true,
        docker = pandas_docker
    }

    # read table with batches information (scatter could be avoided if there was a tail() function)
    Array[Array[String]] mocha_tsv = read_tsv(select_first([mocha_tsv_file]))
    Int n_mocha_batches = length(mocha_tsv)-1
    scatter (idx in range(n_mocha_batches)) { Array[String] mocha_tsv_rows = mocha_tsv[(idx+1)] }
    Map[String, Array[String]] mocha_tbl = as_map(zip(mocha_tsv[0], transpose(mocha_tsv_rows)))
    # check if path is in mocha table (see https://github.com/openwdl/wdl/issues/305)
    Boolean is_path_in_mocha_tbl = length(collect_by_key(zip(flatten([keys(mocha_tbl),["path"]]),range(length(keys(mocha_tbl))+1)))["path"])>1

    # compute data paths for each batch
    scatter (idx in range(n_mocha_batches)) {
      String mocha_data_paths_with_sep = (if defined(mocha_data_path) then sub(select_first([mocha_data_path]), "/$", "") + "/" else "") +
                                         (if is_path_in_mocha_tbl then sub(mocha_tbl["path"][idx], "/$", "") + "/" else "")
    }

    scatter (idx in range(n_mocha_batches)) {
      call vcf_scatter as pgt_scatter {
        input:
          vcf_file = mocha_data_paths_with_sep[idx] + mocha_tbl["pgt_vcf"][idx],
          intervals_bed = ref_scatter_step1.intervals_bed,
          keep_samples_file = prune_file.keep,
          remove_samples_file = remove_samples_file,
          docker = docker_repository_with_sep + bcftools_docker
      }
    }

    Array[Array[File]] interval_slices = transpose(pgt_scatter.vcf_files)
    scatter (idx in range(length(interval_slices))) {
      if (length(interval_slices[idx])>1) {
        call vcf_merge as pgt_merge {
          input:
            vcf_files = interval_slices[idx],
            filebase = filebase + "." + idx,
            docker = docker_repository_with_sep + bcftools_docker
        }
      }

      call pgt_prune {
        input:
          vcf_file = select_first([pgt_merge.vcf_file, interval_slices[idx][0]]),
          sample_tsv_file = sample_tsv_file,
          space_character = space_character,
          min_mac = min_mac_step1,
          min_maf = min_maf_step1,
          autosome_ct = if ref.n_x_chr == 23 then None else ref.n_x_chr - 1,
          docker = docker_repository_with_sep + regenie_docker
      }

      Int n_smpls = select_first([pgt_merge.n_smpls, pgt_scatter.n_smpls[0]])

      if (step1) {
        call regenie_step0 {
          input:
            idx = idx,
            n_phenos = length(select_first([prune_file.pheno_names])),
            n_covars = length(flatten(select_all([prune_file.covar_names]))),
            n_smpls = n_smpls,
            n_markers = pgt_prune.n_markers,
            bed_file = pgt_prune.bed_file,
            bim_file = pgt_prune.bim_file,
            fam_file = pgt_prune.fam_file,
            covar_file = prune_file.covar,
            pheno_file = select_first([prune_file.pheno]),
            binary = binary,
            bsize = bsize_step1,
            loocv = loocv,
            autosome_ct = if ref.n_x_chr == 23 then None else ref.n_x_chr - 1,
            regenie_step0_extra_args = regenie_step0_extra_args,
            filebase = filebase,
            docker = docker_repository_with_sep + regenie_docker
        }
      }
    }

    call pgt_concat {
      input:
        bed_files = pgt_prune.bed_file,
        bim_files = pgt_prune.bim_file,
        fam_files = pgt_prune.fam_file,
        filebase = filebase + ".prune",
        docker = basic_bash_docker
    }

    if (pca && n_smpls[0] >= 50) {
      call plink_pca {
        input:
          n_smpls = n_smpls[0],
          n_markers = pgt_concat.n_markers,
          pca_ndim = pca_ndim,
          pca_cpus = pca_cpus,
          ids_files = pgt_scatter.ids_lines,
          bed_file = pgt_concat.bed_file,
          bim_file = pgt_concat.bim_file,
          fam_file = pgt_concat.fam_file,
          exclusion_regions = ref.pca_exclusion_regions,
          autosome_ct = if ref.n_x_chr == 23 then None else ref.n_x_chr - 1,
          filebase = filebase,
          docker = docker_repository_with_sep + regenie_docker
      }
    }

    if (step1) {
      Array[Array[File]] l0_files = transpose(select_all(regenie_step0.l0_files))
      scatter (idx in range(length(l0_files))) {
        call regenie_step1 {
          input:
            pheno_name = select_first([prune_file.pheno_names])[idx],
            n_covars = length(flatten(select_all([prune_file.covar_names]))),
            n_smpls = n_smpls[0],
            n_markers = pgt_concat.n_markers,
            bed_file = pgt_concat.bed_file,
            bim_file = pgt_concat.bim_file,
            fam_file = pgt_concat.fam_file,
            covar_file = prune_file.covar,
            pheno_file = select_first([prune_file.pheno]),
            n_markers_array = pgt_prune.n_markers,
            l0_files = l0_files[idx],
            binary = binary,
            bsize = bsize_step1,
            loocv = loocv,
            autosome_ct = if ref.n_x_chr == 23 then None else ref.n_x_chr - 1,
            regenie_step1_extra_args = regenie_step1_extra_args,
            filebase = filebase + space_character + select_first([prune_file.pheno_names])[idx],
            docker = docker_repository_with_sep + regenie_docker
        }
        String? loco_lines = if regenie_step1.loco_line == "" then None else regenie_step1.loco_line
        File? loco_file = if defined(loco_lines) then regenie_step1.loco_file else None
        if (binary) {
          String? firth_lines = if select_first([regenie_step1.firth_line]) == "" then None else select_first([regenie_step1.firth_line])
          File? firth_file = if defined(firth_lines) then select_first([regenie_step1.firth_file]) else None
        }
      }
      # unnecessary task for compatibility with Terra https://support.terra.bio/hc/en-us/community/posts/360071465631-write-lines-write-map-write-tsv-write-json-fail-when-run-in-a-workflow-rather-than-in-a-task
      call serialize_lines as loco_lst { input: lines = select_all(loco_lines), filename = filebase + "_pred.list", docker = basic_bash_docker }
      if (binary) {
        call serialize_lines as firth_lst { input: lines = select_all(firth_lines), filename = filebase + "_firth.list", docker = basic_bash_docker }
      }
    }
  }

  if (!step1) {
    if (defined(input_loco_lst) && defined(input_loco_path)) {
      scatter (line in read_lines(select_first([input_loco_lst]))) {
        File input_loco_files = sub(select_first([input_loco_path]), "/$", "") + '/' + sub(line, "^.*[ \t]", "")
      }
    }

    if (binary && defined(input_firth_lst) && defined(input_firth_path)) {
      scatter (line in read_lines(select_first([input_firth_lst]))) {
        File input_firth_files = sub(select_first([input_firth_path]), "/$", "") + '/' + sub(line, "^.*[ \t]", "")
      }
    }
  }

  # REGENIE step 2

  if (step2 || cis) {
    call ref_scatter as ref_scatter_step2 {
      input:
        chrs = select_all(chrs),
        lens = select_all(lens),
        genetic_map_file = ref.genetic_map_file,
        max_win_size_cm = max_win_size_cm_step2,
        overlap_size_cm = 0.0,
        genetic_map_order = false,
        docker = pandas_docker
    }

    Array[Array[String]] intervals_tbl = transpose(read_tsv(ref_scatter_step2.intervals_bed))
    # this is a trick to table how many intervals you will use for each chromosome
    Map[String, Array[Int]] chr_map = collect_by_key(zip(intervals_tbl[0], range(length(intervals_tbl[0]))))

    # read table with batches information (scatter could be avoided if there was a tail() function)
    Array[Array[String]] impute_tsv = read_tsv(select_first([impute_tsv_file]))
    Int n_impute_batches = length(impute_tsv)-1
    scatter (idx in range(n_impute_batches)) { Array[String] impute_tsv_rows = impute_tsv[(idx+1)] }
    Map[String, Array[String]] impute_tbl = as_map(zip(impute_tsv[0], transpose(impute_tsv_rows)))
    # check if path is in impute table (see https://github.com/openwdl/wdl/issues/305)
    Boolean is_path_in_impute_tbl = length(collect_by_key(zip(flatten([keys(impute_tbl),["path"]]),range(length(keys(impute_tbl))+1)))["path"])>1

    # compute data paths for each batch
    scatter (idx in range(n_impute_batches)) {
      String impute_data_paths_with_sep = (if defined(impute_data_path) then sub(select_first([impute_data_path]), "/$", "") + "/" else "") +
                                          (if is_path_in_impute_tbl then sub(impute_tbl["path"][idx], "/$", "") + "/" else "")
    }

    scatter (p in cross(range(n_impute_batches), range(length(select_all(chrs))))) {
      File imp_vcf_file = impute_data_paths_with_sep[p.left] + impute_tbl[("chr" + sub(select_all(chrs)[p.right], "^chr", "") + "_imp_vcf")][p.left]
      if (length(chr_map[(select_all(chrs)[p.right])]) > 1) {
        call vcf_scatter {
          input:
            vcf_file = imp_vcf_file,
            intervals_bed = ref_scatter_step2.intervals_bed,
            keep_samples_file = prune_file.keep,
            remove_samples_file = remove_samples_file,
            chr = select_all(chrs)[p.right],
            dosage_field = dosage_field,
            docker = docker_repository_with_sep + bcftools_docker
        }
      }
      Int cross_idx = p.right
      Array[File] scatter_vcf_files = select_first([vcf_scatter.vcf_files, [imp_vcf_file]])
    }

    Map[Int, Array[Array[File]]] idx2vcf_files = collect_by_key(zip(cross_idx, scatter_vcf_files))
    scatter (idx in range(length(select_all(chrs)))) { Array[Array[File]] slices_vcf_files = transpose(idx2vcf_files[idx]) }
    Array[Array[File]] matrix_vcf_files = flatten(slices_vcf_files)

    if (step2) {
      # generate list of expected output association files
      scatter (line in read_lines(select_first([firth_lst.file, loco_lst.file, input_firth_lst, input_loco_lst]))) {
        String regenie_suffix = sub(line, " .*$", "") + (if defined(pop) then "." + select_first([pop]) else "") + ".regenie.gz" # https://github.com/broadinstitute/cromwell/issues/5549
      }
    }

    # merging has to happen at the VCF level as plink2 does not currently merge pgen files
    # https://www.cog-genomics.org/plink/2.0/data#pmerge
    scatter (idx in range(length(matrix_vcf_files))) {
      if (length(matrix_vcf_files[idx])>1 || defined(min_mac_step2)) {
        call vcf_merge {
          input:
            vcf_files = matrix_vcf_files[idx],
            min_mac = min_mac_step2,
            filebase = filebase + "." + idx,
            docker = docker_repository_with_sep + bcftools_docker
        }
      }

      call vcf2pgen {
        input:
          vcf_file = select_first([vcf_merge.vcf_file, matrix_vcf_files[idx][0]]),
          sample_tsv_file = sample_tsv_file,
          dosage_field = dosage_field,
          space_character = space_character,
          autosome_ct = if ref.n_x_chr == 23 then None else ref.n_x_chr - 1,
          par_bp1 = ref.par_bp1,
          par_bp2 = ref.par_bp2,
          docker = docker_repository_with_sep + regenie_docker
      }

      if (step2) {
        call regenie_step2 {
          input:
            chr = intervals_tbl[0][idx],
            n_phenos = length(select_first([prune_file.pheno_names])),
            n_covars = length(flatten(select_all([prune_file.covar_names]))),
            n_smpls = vcf2pgen.n_smpls,
            n_markers = vcf2pgen.n_markers,
            fasta_fai = ref.fasta_fai,
            pgen_file = vcf2pgen.pgen_file,
            pvar_file = vcf2pgen.pvar_file,
            psam_file = vcf2pgen.psam_file,
            covar_file = prune_file.covar,
            pheno_file = select_first([prune_file.pheno]),
            pop = pop,
            regenie_suffix = select_first([regenie_suffix]), # https://github.com/broadinstitute/cromwell/issues/5549
            binary = binary,
            bsize = bsize_step2,
            min_info = min_info_step2,
            autosome_ct = if ref.n_x_chr == 23 then None else ref.n_x_chr - 1,
            regenie_step2_extra_args = regenie_step2_extra_args,
            loco_lst = select_first([loco_lst.file, input_loco_lst]),
            loco_files = select_first([loco_files, input_loco_files]),
            firth_lst = if binary then select_first([firth_lst.file, input_firth_lst]) else None,
            firth_files = if binary then select_first([firth_files, input_firth_files]) else None,
            docker = docker_repository_with_sep + regenie_docker
        }
      }
    }

    if (step2) {
      call vcf_concat {
        input:
          vcf_files = select_all(regenie_step2.vcf_file),
          ref_fasta = if defined(ref.gff3_file) then ref.fasta else None,
          fasta_fai = if defined(ref.gff3_file) then ref.fasta_fai else None,
          gff3_file = ref.gff3_file,
          rsid_vcf_file = ref.rsid_vcf_file,
          rsid_vcf_idx = ref.rsid_vcf_idx,
          filebase = filebase + (if length(select_first([prune_file.pheno_names])) == 1 then "." + select_first([prune_file.pheno_names])[0] else "") + (if defined(pop) then "." + select_first([pop]) else "") + ".gwas",
          docker = docker_repository_with_sep + bcftools_docker
      }
      Array[Array[File]] regenie_matrix_files = transpose(select_all(regenie_step2.regenie_files))
      scatter (idx in range(length(regenie_matrix_files))) {
        call assoc_concat as regenie_concat {
          input:
            assoc_files = regenie_matrix_files[idx],
            n_x_chr = ref.n_x_chr,
            filebase = filebase + "." + select_first([regenie_suffix])[idx],
            docker = docker_repository_with_sep + bcftools_docker
        }
        if (plot && regenie_concat.has_data) {
            call assoc_plot as regenie_plot {
              input:
                assoc_file = regenie_concat.file,
                genome = if ref_name == "GRCh38" || ref_name == "GRCh37" then ref_name else None,
                cyto_file = ref.cyto_file,
                autosome_ct = if ref.n_x_chr == 23 then None else ref.n_x_chr - 1,
                filebase = basename(filebase + "." + select_first([regenie_suffix])[idx], ".gz"),
                docker = docker_repository_with_sep + r_mocha_docker
          }
        }
      }
    }

    if (cis) {
      Array[String]+ lines = if defined(loco_lst.file) || defined(input_loco_lst) then read_lines(select_first([loco_lst.file, input_loco_lst])) else select_first([prune_file.pheno_names])
      scatter (idx in range(length(lines))) {
        String plink_pheno_names = if defined(loco_lst.file) || defined(input_loco_lst) then sub(lines[idx], " .*$", "") else lines[idx]
        String plink_pheno_chrs = sub(sub(sub(plink_pheno_names, "_.*$", ""), "[pq]*$", ""), "Y", "X")
        # check if plink_pheno_chrs is present in the chr2idx to know whether the cis association should be run (see https://github.com/openwdl/wdl/issues/305)
        Int? cis_idx = if length(collect_by_key(zip(flatten([keys(chr2idx),[plink_pheno_chrs]]),range(length(keys(chr2idx))+1)))[plink_pheno_chrs])>1 then idx else None
      }
      # this map, given a chromosome (1, 2, ..., X), returns the indexes of the intervals for that chromosomes
      scatter (chr in intervals_tbl[0]) { String chr_string = sub(chr, "^chr", "") }
      Map[String, Array[Int]] chr2idx = collect_by_key(zip(chr_string, range(length(intervals_tbl[0]))))
      # the following code checks that the chromosome name is in the list of available chromosomes
      scatter (idx in select_all(cis_idx)) {
        Array[Pair[Int, Int]] pheno_interval_pairs = cross([idx], chr2idx[(plink_pheno_chrs[idx])])
      }
      # maybe I should test whether the interval falls under the event or not
      String x_chr_num = ref.n_x_chr
      scatter (p in flatten(pheno_interval_pairs)) {
        Int pheno_idx = p.left
        call plink_glm {
          input:
            chr_num = sub(plink_pheno_chrs[p.left], "X", x_chr_num),
            pheno_name = plink_pheno_names[p.left],
            n_phenos = length(select_first([prune_file.pheno_names])),
            n_covars = length(flatten(select_all([prune_file.covar_names]))),
            n_smpls = vcf2pgen.n_smpls[p.right],
            n_markers = vcf2pgen.n_markers[p.right],
            pgen_file = vcf2pgen.pgen_file[p.right],
            pvar_file = vcf2pgen.pvar_file[p.right],
            psam_file = vcf2pgen.psam_file[p.right],
            loco_file =  if defined(loco_files) || defined(input_loco_files) then select_first([loco_files, input_loco_files])[p.left] else None,
            covar_file = prune_file.covar,
            pheno_file = select_first([prune_file.pheno]),
            binary = binary,
            max_vif = max_vif,
            max_corr = max_corr,
            plink_extra_args = plink_extra_args,
            docker = docker_repository_with_sep + regenie_docker
        }
      }

      Map[Int, Array[File]] idx2assoc_files = collect_by_key(zip(pheno_idx, plink_glm.assoc_file))
      scatter (idx in select_all(cis_idx)) {
        call assoc_concat as plink_concat {
          input:
            assoc_files = idx2assoc_files[idx],
            n_x_chr = ref.n_x_chr,
            zst = true,
            filebase = filebase + "." + plink_pheno_names[idx] + ".glm." + (if binary then "logistic.hybrid" else "linear") + ".gz",
            docker = docker_repository_with_sep + regenie_docker
        }
        if (plot && plink_concat.has_data) {
          call assoc_plot as plink_plot {
            input:
              assoc_file = plink_concat.file,
              genome = if ref_name == "GRCh38" || ref_name == "GRCh37" then ref_name else None,
              cyto_file = ref.cyto_file,
              autosome_ct = if ref.n_x_chr == 23 then None else ref.n_x_chr - 1,
              min_af = cis_plot_min_af,
              filebase = filebase + "." + plink_pheno_names[idx] + ".glm." + (if binary then "logistic.hybrid" else "linear"),
              docker = docker_repository_with_sep + r_mocha_docker
          }
        }
      }
    }
  }

  output {
    File? bed_file = if pgt_output then pgt_concat.bed_file else None
    File? bim_file = if pgt_output then pgt_concat.bim_file else None
    File? fam_file = if pgt_output then pgt_concat.fam_file else None
    File? eigenvec_file = plink_pca.eigenvec_file
    File? eigenval_file = plink_pca.eigenval_file
    File? pcs_tsv_file = plink_pca.pcs_tsv_file
    File? loco_lst_file = loco_lst.file
    Array[File]? loco_files = if defined(loco_file) then select_all(select_first([loco_file])) else None
    File? firth_lst_file = firth_lst.file
    Array[File]? firth_files = if defined(firth_file) then select_all(select_first([firth_file])) else None
    File? gwas_vcf_file = vcf_concat.vcf_file
    File? gwas_vcf_idx = vcf_concat.vcf_idx
    Array[File]? regenie_files = regenie_concat.file
    Array[File]? regenie_indexes = regenie_concat.index
    Array[File]? regenie_png_files = if step2 && plot then select_all(select_first([regenie_plot.png_file])) else None
    Array[File]? plink_files = plink_concat.file
    Array[File]? plink_indexes = plink_concat.index
    Array[File]? plink_png_files = if cis && plot then select_all(select_first([plink_plot.png_file])) else None
    Array[File]? regenie_step0_logs = if step1 then select_all(select_first([regenie_step0.log_file])) else None
    Array[File]? regenie_step1_logs = if step1 then select_all(select_first([regenie_step1.log_file])) else None
    Array[File]? regenie_step2_logs = if step2 then select_all(select_first([regenie_step2.log_file])) else None
    Array[File]? plink_logs = plink_glm.log_file
  }

  meta {
    author: "Giulio Genovese"
    email: "giulio.genovese@gmail.com"
    description: "See the [MoChA](https://github.com/freeseek/mocha) website for more information"
  }
}

task get_n {
  input {
    File file

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mv "~{file}" .
    grep -v ^# "~{basename(file)}" | wc -l
    rm "~{basename(file)}"
  >>>

  output {
    Int n = read_int(stdout())
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

# use of !(a!=b) due to bug Cromwell team will not fix: https://github.com/broadinstitute/cromwell/issues/5602
task prune_file {
  input {
    String? sex_specific
    File sample_tsv_file
    File? keep_samples_file
    File? remove_samples_file
    File? covar_tsv_file
    File pheno_tsv_file
    String space_character
    Boolean binary
    Int min_case_count
    Int min_sex_count
    String filebase

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    echo "~{sep("\n", select_all([sample_tsv_file, keep_samples_file, remove_samples_file, covar_tsv_file, pheno_tsv_file]))}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    awk -F"\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
      NR>1 {sex=substr($(f["computed_gender"]),1,1); if (toupper(sex)=="M" || sex==1) printf "%s\t1\n",$(f["sample_id"])}' \
      "~{basename(sample_tsv_file)}" > "~{filebase}.male"
    awk -F"\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
      NR>1 {sex=substr($(f["computed_gender"]),1,1); if (toupper(sex)=="F" || sex==2) printf "%s\t2\n",$(f["sample_id"])}' \
      "~{basename(sample_tsv_file)}" > "~{filebase}.female"
    ~{if defined(sex_specific) || defined(keep_samples_file) then "cut -f1 \"" + filebase + "." +
      (if defined(sex_specific) then select_first([sex_specific]) else "male\" \"" + filebase + ".female") + "\"" +
      (if defined(keep_samples_file) then " | \\\n  awk -F\"\\t\" 'NR==FNR {x[$1]++} NR>FNR && $1 in x' \"" + basename(select_first([keep_samples_file])) + "\" -" else "") +
      (if defined(remove_samples_file) then " | \\\n  awk -F\"\\t\" 'NR==FNR {x[$1]++} NR>FNR && !($1 in x)' \"" + basename(select_first([remove_samples_file])) + "\" -" else "") +
      " > \"" + filebase + ".keep.lines\"\n"
      else ""}cat "~{filebase + "." + if defined(sex_specific) then select_first([sex_specific]) else "male\" \"" + filebase + ".female"}" | \
      ~{if defined(remove_samples_file) then
        "awk -F\"\\t\" 'NR==FNR {x[$1]++} NR>FNR && !($1 in x)' \"" + basename(select_first([remove_samples_file])) + "\" - | \\\n  "
        else "" + if defined(covar_tsv_file) then
        "awk -F\"\\t\" 'NR==FNR {x[$1]++} NR>FNR && $1 in x' \"" + basename(select_first([covar_tsv_file])) + "\" - | \\\n  "
        else ""}awk -F"\t" 'NR==FNR {x[$1]++} NR>FNR && (FNR==1 || $1 in x)' - "~{basename(pheno_tsv_file)}" > "~{filebase}.tmp"
    cat "~{filebase}.male" "~{filebase}.female" | \
      awk -F"\t" 'NR==FNR {sex[$1]=$2} NR>FNR && FNR==1 {for (i=2; i<=NF; i++) pheno[i] = $i}
      NR>FNR && FNR>1 {for (i=2; i<=NF; i++) {if ($i==0) ctrls[i]++; if ($i==1) cases[i]++
      if (sex[$1]==1 && $i!="NA") males[i]++; if (sex[$1]==2 && $i!="NA") females[i]++}}
      END {for (i in pheno) printf "%s\t%d\t%d\t%d\t%d\n",pheno[i],ctrls[i],cases[i],males[i],females[i]}' \
      - "~{filebase}.tmp" > "~{filebase}.cnt"
    awk -F"\t" 'NR==FNR ~{if binary then "&& $2>=" + min_case_count + " && $3>=" + min_case_count else ""} && $~{
      if defined(sex_specific) && !(select_first([sex_specific]) != "male") then "4"
      else if defined(sex_specific) && !(select_first([sex_specific]) != "female") then "5"
      else "4>=" + min_sex_count + " && $5"}>=~{min_sex_count} {keep[$1]++}
      NR>FNR {if (FNR==1) {for (i=2; i<=NF; i++) if ($i in keep) col[j++]=i; printf "FID\tIID"}
      else {gsub(" ","~{space_character}",$1); printf "0\t%s",$1} for (i=0; i<j; i++) printf "\t%s",$col[i]; printf "\n"}' \
      "~{filebase}.cnt" "~{filebase}.tmp" > "~{filebase}.phe"
    ~{if defined(covar_tsv_file) then "cat \"" + filebase + "." +
      (if defined(sex_specific) then select_first([sex_specific]) else "male\" \"" + filebase + ".female") + "\" | \\\n" +
      "  awk -F\"\\t\" 'NR==FNR {sex[$1]=$2} NR>FNR && (FNR==1 || $1 in sex) {if (FNR==1) printf \"FID\\tIID" + (if defined(sex_specific) then "" else "\\tsex") + "\"\n" +
      "  else {gsub(\" \",\"" + space_character + "\",$1); printf \"0\\t%s" + (if defined(sex_specific) then "" else "\\t%s") + "\",$1" + (if defined(sex_specific) then "" else ",sex[$1]") + "}\n" +
      "  for (i=2; i<=NF; i++) printf \"\\t%s\",$i; printf \"\\n\"}' \\\n" +
      "  - \"" + basename(select_first([covar_tsv_file])) + "\" > \"" + filebase + ".cov\"\n"
      else if !defined(sex_specific) then
      "cat \"" + filebase + ".male\" \"" + filebase + ".female\" | \\\n" +
      "  awk -F\"\\t\" '{if (NR==1) print \"FID\\tIID\\tsex\"; else printf \"0\\t%s\\t%s\\n\",$1,$2}' > \"" + filebase + ".cov\"\n"
      else ""}head -n1 "~{filebase}.phe" | cut -f3- | tr '\t' '\n'
    ~{if defined(covar_tsv_file) || !defined(sex_specific) then
      "head -n1 \"" + filebase + ".cov\" | cut -f3- | tr '\\t' '\\n' > \"" + filebase + ".cov.lines\"\n"
    else ""}rm "~{filebase}.male" "~{filebase}.female" "~{filebase}.tmp" "~{filebase}.cnt"
    echo "~{sep("\n", select_all([sample_tsv_file, keep_samples_file, remove_samples_file, covar_tsv_file, pheno_tsv_file]))}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    Array[String] pheno_names = read_lines(stdout())
    Array[String]? covar_names = if defined(covar_tsv_file) || !defined(sex_specific) then read_lines(filebase + ".cov.lines") else None
    File? keep = if defined(sex_specific) || defined(keep_samples_file) then filebase + ".keep.lines" else None
    File pheno = filebase + ".phe"
    File? covar = if defined(covar_tsv_file) || !defined(sex_specific) then filebase + ".cov" else None
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

task ref_scatter {
  input {
    Array[String]+ chrs
    Array[String]+ lens
    File genetic_map_file
    Float max_win_size_cm
    Float overlap_size_cm
    Boolean genetic_map_order

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mv "~{genetic_map_file}" .
    chrs=~{write_lines(chrs)}
    lens=~{write_lines(lens)}
    paste -d $'\t' $chrs $lens > chr2len.tsv
    python3 <<CODE
    import sys, pandas as pd, numpy as np
    chr2len = {}
    with open('chr2len.tsv') as f:
      for line in f:
        (key, val) = line.split('\t')
        chr2len[key] = int(val)
    df_map = pd.read_csv('~{basename(genetic_map_file)}', delim_whitespace = True, header = 0, names = ['CHR', 'POS' ,'RATE', 'CM'])
    df_out = {}
    for chr, df_group in df_map.groupby('CHR'):
      fai_chr = str(chr) if str(chr) in chr2len else 'chr' + str(chr) if 'chr' + str(chr) in chr2len else 'X' if 'X' in chr2len else 'chrX' if 'chrX' in chr2len else None
      if fai_chr:
        chr_cm_len = max(df_group['CM'])
        n_win = np.ceil((chr_cm_len - ~{overlap_size_cm})/(~{max_win_size_cm} - ~{overlap_size_cm}))
        win_size = (chr_cm_len - ~{overlap_size_cm}) / n_win + ~{overlap_size_cm}
        cm_begs = (win_size - ~{overlap_size_cm}) * np.arange(1, n_win)
        cm_ends = (win_size - ~{overlap_size_cm}) * np.arange(1, n_win) + ~{overlap_size_cm}
        pos_begs = np.concatenate(([0], 0 + np.interp(cm_begs, df_group['CM'], df_group['POS'], period = np.inf).astype(int)))
        pos_ends = np.concatenate((np.interp(cm_ends, df_group['CM'], df_group['POS'], period = np.inf).astype(int), [chr2len[fai_chr]]))
        df_out[fai_chr] = pd.DataFrame.from_dict({'CHR': fai_chr, 'BEG': pos_begs, 'END': pos_ends})
    df = pd.concat(~{if genetic_map_order then "df_out" else "[df_out[fai_chr] for fai_chr in chr2len.keys()]"})
    df[['CHR', 'BEG', 'END']].to_csv('ref_scatter.bed', sep='\t', header = False, index = False)
    CODE
    rm chr2len.tsv
    rm "~{basename(genetic_map_file)}"
  >>>

  output {
    File intervals_bed = "ref_scatter.bed"
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

# the command requires BCFtools 1.15 due to bug https://github.com/samtools/bcftools/issues/1631
task vcf_scatter {
  input {
    File vcf_file
    File intervals_bed # zero-based intervals
    File? keep_samples_file
    File? remove_samples_file
    Int clevel = 2
    String? chr
    String? dosage_field

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 3.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])
  String filebase = basename(basename(vcf_file, ".bcf"), ".vcf.gz")

  command <<<
    set -euo pipefail
    echo "~{sep("\n", select_all([vcf_file, intervals_bed, keep_samples_file, remove_samples_file]))}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    ~{if defined(chr) then
      "mv \"" + basename(intervals_bed) + "\" \"" + basename(intervals_bed, ".bed") + ".all.bed\"\n" +
      "awk -v chr=\"" + chr + "\" '$1==chr' \"" + basename(intervals_bed, ".bed") + ".all.bed\" > \"" + basename(intervals_bed) + "\"\n" +
      "rm \"" + basename(intervals_bed, ".bed") + ".all.bed\""
      else ""}
    bcftools query --force-samples --list-samples ~{if defined(keep_samples_file) then
      "--samples-file \"" + basename(select_first([keep_samples_file])) + "\" "
    else if defined (remove_samples_file) then
      "--samples-file \"^" + basename(select_first([remove_samples_file])) + "\" "
    else ""}"~{basename(vcf_file)}" > "~{filebase}.ids.lines"
    cat "~{filebase}.ids.lines" | wc -l > n_smpls.int
    awk -F"\t" '{print $1":"1+$2"-"$3"\t"NR-1}' "~{basename(intervals_bed)}" > regions.lines
    ~{if defined(keep_samples_file) then
      "bcftools view --no-version -Ou --samples-file \"" + basename(select_first([keep_samples_file])) + "\" --force-samples \"" + basename(vcf_file) + "\" |\n  "
    else if defined (remove_samples_file) then
      "bcftools view --no-version -Ou --samples-file \"^" + basename(select_first([remove_samples_file])) + "\" --force-samples \"" + basename(vcf_file) + "\" |\n  "
    else ""}bcftools annotate \
      --no-version \
      --output-type u \
      --remove ID,QUAL,FILTER,INFO,^FMT/GT~{if defined(dosage_field) then ",^FMT/DS" else ""} \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      ~{if defined(keep_samples_file) || defined(remove_samples_file) then "" else "\"" + basename(vcf_file) + "\" "} | \
    bcftools norm \
      --no-version \
      --output-type u \
      --rm-dup exact \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} | \
    bcftools +scatter \
      --no-version \
      --output-type b~{clevel} \
      --output vcfs \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      --scatter-file regions.lines \
      --prefix "~{filebase}."
    cut -f2 regions.lines | sed 's/^/vcfs\/~{filebase}./;s/$/.bcf/'
    echo "~{sep("\n", select_all([vcf_file, intervals_bed, keep_samples_file, remove_samples_file]))}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
    rm regions.lines
  >>>

  output {
    Int n_smpls = read_int("n_smpls.int")
    File ids_lines = filebase + ".ids.lines"
    Directory vcfs = "vcfs"
    Array[File] vcf_files = read_lines(stdout())
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

task vcf_merge {
  input {
    Array[File]+ vcf_files
    Int? min_mac
    Int clevel = 2
    String filebase

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10 + 2.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])

  command <<<
    set -euo pipefail
    ~{if length(vcf_files) > 1 then "vcf_files=" else ""}~{if length(vcf_files) > 1 then write_lines(vcf_files) else ""}
    ~{if length(vcf_files) > 1 then "cat $vcf_files | tr '\\n' '\\0' | xargs -0 mv -t .\n" +
      "sed -i 's/^.*\\///' $vcf_files\n" +
      "bcftools merge \\\n" +
      "  --no-version \\\n" +
      "  --output-type " + (if defined(min_mac) then "u" else "b" + clevel + " \\\n" +
      "  --output \"" + filebase + ".bcf\"") + " \\\n" +
      "  --file-list $vcf_files \\\n" +
      "  --merge none \\\n" +
      "  --no-index \\\n" +
      (if cpu > 1 then "  --threads " + (cpu - 1) else "") +
      (if defined(min_mac) then " | \\\n" +
      "bcftools view \\\n" +
      "  --no-version \\\n" +
      "  --output-type b" + clevel + " \\\n" +
      "  --output \"" + filebase + ".bcf\" \\\n" +
      "  --min-ac " + min_mac + ":nonmajor" else "")
      else if defined(min_mac) then
      "bcftools view \\\n" +
      "  --no-version \\\n" +
      "  --output-type b" + clevel + " \\\n" +
      "  --output \"" + filebase + ".bcf\" \\\n" +
      "  --min-ac " + min_mac + ":nonmajor \\\n" +
      "  \"" + vcf_files[0] + "\""
      else "mv \"" + vcf_files[0] + "\" \"" + filebase + ".bcf\""}
    bcftools query --list-samples "~{filebase}.bcf" | wc -l
    ~{if length(vcf_files) > 1 then "cat $vcf_files | tr '\\n' '\\0' | xargs -0 rm" else ""}
  >>>

  output {
    File vcf_file = filebase + ".bcf"
    Int n_smpls = read_int(stdout())
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

# this command needs PLINK 1.9 as conversion from VCF using PLINK 2.0 is inefficient: https://groups.google.com/g/plink2-users/c/hsByNOklyA0
# the U sex needs to be encoded as 0 as this is the only accepted value for PLINK 1.9: https://groups.google.com/g/plink2-users/c/z7YJYa677NQ
# use of !(a!=b) due to bug Cromwell team will not fix: https://github.com/broadinstitute/cromwell/issues/5602
# the command requires BCFtools 1.14 due to bug https://github.com/samtools/bcftools/issues/1528
task pgt_prune {
  input {
    File vcf_file
    File sample_tsv_file
    String space_character
    Int min_mac
    Float min_maf
    Int? autosome_ct

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  String filebase = basename(basename(vcf_file, ".bcf"), ".vcf.gz")
  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 3.0 * vcf_size)])

  command <<<
    set -euo pipefail
    mv "~{vcf_file}" .
    mv "~{sample_tsv_file}" .
    awk 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
      NR>1 {id=$(f["sample_id"]); gsub(" ","~{space_character}",id);
      print 0,id,toupper(substr($(f["computed_gender"]),1,1))}' "~{basename(sample_tsv_file)}" | \
      sed 's/U$/0/;s/K$/1/' > "~{filebase}.sex"
    rm "~{basename(sample_tsv_file)}"
    bcftools +fill-tags --no-version -Ou --include 'sum(AC)>=~{min_mac} && AN-sum(AC)>=~{min_mac} && MAF>=~{min_maf}' "~{basename(vcf_file)}" -- --tags AC,AN,MAF | \
    bcftools annotate --no-version -Ob0 --set-id "%VKX" --remove FILTER,INFO,^FMT/GT | \
    plink1.9 \
      --threads ~{cpu} \
      --memory ~{round(1024 * memory - 512)} \
      --bcf /dev/stdin \
      --update-sex "~{filebase}.sex" \
      --keep-allele-order \
      --vcf-idspace-to ~{space_character} \
      --const-fid \
      --allow-extra-chr 0 \
      ~{if defined(autosome_ct) then "--chr-set " + autosome_ct else ""} \
      --make-bed \
      --out "~{filebase}" \
      1>&2
    rm "~{basename(vcf_file)}" "~{filebase}.sex" "~{filebase}.nosex"
    plink1.9 \
      --threads ~{cpu} \
      --memory ~{round(1024 * memory - 512)} \
      --bfile "~{filebase}" \
      --keep-allele-order \
      --indep 50 5 2 \
      ~{if defined(autosome_ct) then "--chr-set " + autosome_ct else ""} \
      --out "~{filebase}" \
      1>&2
    rm "~{filebase}.prune.out"
    cat "~{filebase}.prune.in" | wc -l
    plink1.9 \
      --memory ~{round(1024 * memory - 512)} \
      --threads ~{cpu} \
      --bfile "~{filebase}" \
      --keep-allele-order \
      --extract "~{filebase}.prune.in" \
      ~{if defined(autosome_ct) then "--chr-set " + autosome_ct else ""} \
      --make-bed \
      --out "~{filebase}.prune" \
      1>&2
    rm "~{filebase}.bed" "~{filebase}.bim" "~{filebase}.fam" "~{filebase}.prune.in"
  >>>

  output {
    Int n_markers = read_int(stdout())
    File bed_file = filebase + ".prune.bed"
    File bim_file = filebase + ".prune.bim"
    File fam_file = filebase + ".prune.fam"
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

task pgt_concat {
  input {
    Array[File]+ bed_files
    Array[File]+ bim_files
    Array[File]+ fam_files
    String filebase

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 0
    Int maxRetries = 0
  }

  Float bed_size = size(bed_files, "GiB")
  Float bim_size = size(bim_files, "GiB")
  Float fam_size = size(fam_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * (bed_size + bim_size + fam_size))])

  command <<<
    set -euo pipefail
    bed_files=~{write_lines(bed_files)}
    bim_files=~{write_lines(bim_files)}
    fam_files=~{write_lines(fam_files)}
    cat $bed_files $bim_files $fam_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $bed_files $bim_files $fam_files
    (echo -en "\x6C\x1B\x01"; cat $bed_files | tr '\n' '\0' | xargs -0 tail -qc+4) > "~{filebase}.bed"
    cat $bim_files | tr '\n' '\0' | xargs -0 cat > "~{filebase}.bim"
    head -n1 $fam_files | tr '\n' '\0' | xargs -0 cat > "~{filebase}.fam"
    cat "~{filebase}.bim" | wc -l
    cat $bed_files $bim_files $fam_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    Int n_markers = read_int(stdout())
    File bed_file = filebase + ".bed"
    File bim_file = filebase + ".bim"
    File fam_file = filebase + ".fam"
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

task plink_pca {
  input {
    Int n_smpls
    Int n_markers
    Int pca_ndim = 20
    Int pca_cpus = 2
    Array[File]+ ids_files
    File bed_file
    File bim_file
    File fam_file
    String? exclusion_regions
    Int? autosome_ct
    String filebase

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 0
    Int maxRetries = 0
  }

  Float bed_size = size(bed_file, "GiB")
  Float bim_size = size(bim_file, "GiB")
  Float fam_size = size(fam_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * (bed_size + bim_size + fam_size))])
  # https://groups.google.com/g/plink2-users/c/qGiWkqhuvcY/m/kpeS3FfQBAAJ
  Float memory = select_first([memory_override, 3.5 + (16.0 * pca_ndim * (pca_ndim + 1) * (n_markers + if n_markers > n_smpls then n_markers else n_smpls) + (16.0 * pca_ndim * (pca_ndim + 1) + 5760.0) * n_smpls * pca_cpus) / 1024 / 1024 / 1024])
  Int cpu = select_first([cpu_override, if 2 * ceil(memory / 13) > pca_cpus then 2 * ceil(memory / 13) else pca_cpus]) # always require at least two CPUs

  command <<<
    set -euo pipefail
    ids_files=~{write_lines(ids_files)}
    echo "~{sep("\n", [bed_file, bim_file, fam_file])}" | \
      cat $ids_files - | tr '\n' '\0' | xargs -0 mv -t .
    cat $ids_files | sed 's/^.*\///' | tr '\n' '\0' | xargs -0 cat > ids.lines
    ~{if defined(exclusion_regions) then "echo \"" + exclusion_regions + "\" | \\\n" +
      "tr ',[:\\-]' '\\n ' | awk '{print $0,\"r\"NR}' > exclusion_regions.txt"
    else ""}
    plink2 \
      --threads ~{cpu} \
      --memory ~{round(1024 * memory - 512)} \
      --bed "~{basename(bed_file)}" \
      --bim "~{basename(bim_file)}" \
      --fam "~{basename(fam_file)}" \
      ~{if defined(exclusion_regions) then "--exclude range exclusion_regions.txt" else ""} \
      --pca ~{pca_ndim} approx \
      ~{if defined(autosome_ct) then "--chr-set " + autosome_ct else ""} \
      --out "~{filebase}"
    ~{if defined(exclusion_regions) then "rm exclusion_regions.txt" else ""}
    (head -n1 "~{filebase}.eigenvec" | cut -f2- | sed 's/^IID/sample_id/;s/PC/pc/g';
    awk -F"\t" -v OFS="\t" 'NR==FNR {e[NR]=$1} NR>FNR && FNR>1 {for (i=3; i<=NF; i++) $i*=e[i-2]; print}' \
      "~{filebase}.eigenval" "~{filebase}.eigenvec" | cut -f3- | paste -d $'\t' ids.lines -) > "~{filebase}.pcs.tsv"
    rm ids.lines
    echo "~{sep("\n", [bed_file, bim_file, fam_file])}" | \
      cat $ids_files - | sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
    set -euo pipefail
  >>>

  output {
    File eigenvec_file = filebase + ".eigenvec"
    File eigenval_file = filebase + ".eigenval"
    File pcs_tsv_file = filebase + ".pcs.tsv"
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

# see https://github.com/rgcgithub/regenie/wiki/Further-parallelization-for-level-0-models-in-Step-1
task regenie_step0 {
  input {
    Int idx
    Int n_phenos
    Int n_covars
    Int n_smpls
    Int n_markers
    Int n_ridge_l0 = 5
    File bed_file
    File bim_file
    File fam_file
    File? covar_file
    File pheno_file
    Boolean binary
    Int bsize
    Boolean loocv
    Int? autosome_ct
    String? regenie_step0_extra_args
    String filebase

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0

    Float mult = 2.0 # this should not be necessary but most likely there is be a bug in REGENIE
  }

  Float bed_size = size(bed_file, "GiB")
  Float bim_size = size(bim_file, "GiB")
  Float fam_size = size(fam_file, "GiB")
  # see print_usage_info() in https://github.com/rgcgithub/regenie/blob/master/src/Regenie.cpp
  Int disk_size = select_first([disk_size_override, ceil(10.0 + bed_size + bim_size + fam_size + 8.0 * n_phenos * ceil(n_markers / bsize) * n_ridge_l0 * n_smpls / 1024 / 1024 / 1024)])
  Float memory = select_first([memory_override, 3.5 + mult * 8.0 * (bsize + n_phenos * (4.0 + n_ridge_l0) + n_covars) * n_smpls / 1024 / 1024 / 1024])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])
  Int n_bins = ceil(1.0 * n_markers / bsize)

  command <<<
    set -euo pipefail
    echo "~{sep("\n", select_all([bed_file, bim_file, fam_file, covar_file, pheno_file]))}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    echo -e "~{n_markers} ~{bsize}\n~{filebase}_job~{idx+1} ~{n_bins} ~{n_markers}" > "~{filebase}_~{idx+1}.master"
    cut -f2 "~{basename(bim_file)}" > "~{filebase}_job~{idx+1}.snplist"
    regenie \
      --step 1 \
      --bed "~{basename(bed_file, ".bed")}" \
      --phenoFile "~{basename(pheno_file)}" \
      ~{if defined(covar_file) then "--covarFile \"" + basename(select_first([covar_file])) + "\"" else ""} \
      --gz \
      ~{if binary then "--bt" else ""} \
      --bsize ~{bsize} \
      ~{if loocv then "--loocv" else ""} \
      --run-l0 "~{filebase}_~{idx+1}.master",1 \
      ~{if defined(autosome_ct) then "--nauto " + autosome_ct else ""} \
      ~{if defined(regenie_step0_extra_args) then regenie_step0_extra_args else ""} \
      --threads ~{cpu} \
      --out "~{filebase}" \
      1>&2
    rm "~{filebase}_~{idx+1}.master" "~{filebase}_job~{idx+1}.snplist"
    mkdir logs
    mv "~{filebase}.log" "logs/~{filebase}.~{idx+1}.step0.log"
    mkdir l0
    seq 1 ~{n_phenos} | sed 's/^/~{filebase}_job~{idx+1}_l0_Y/' | tr '\n' '\0' | xargs -0 mv -t l0
    seq 1 ~{n_phenos} | sed 's/^/l0\/~{filebase}_job~{idx+1}_l0_Y/'
    echo "~{sep("\n", select_all([bed_file, bim_file, fam_file, covar_file, pheno_file]))}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File log_file = "logs/" + filebase + "." + (idx + 1) + ".step0.log"
    Directory l0_dir = "l0"
    Array[File] l0_files = read_lines(stdout())
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

# see https://github.com/rgcgithub/regenie/wiki/Further-parallelization-for-level-0-models-in-Step-1
task regenie_step1 {
  input {
    String pheno_name
    Int n_covars
    Int n_smpls
    Int n_markers
    File bed_file
    File bim_file
    File fam_file
    File? covar_file
    File pheno_file
    Array[Int]+ n_markers_array
    Array[File]+ l0_files
    Boolean binary
    Int bsize
    Boolean loocv
    Int? autosome_ct
    String? regenie_step1_extra_args
    String filebase

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 0
    Int maxRetries = 0
  }

  Float bed_size = size(bed_file, "GiB")
  Float bim_size = size(bim_file, "GiB")
  Float fam_size = size(fam_file, "GiB")
  Float l0_size = size(l0_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + bed_size + bim_size + fam_size + l0_size)])
  Float memory = select_first([memory_override, 3.5 + bed_size + bim_size + fam_size + l0_size])
  Int cpu = select_first([cpu_override, 2 * ceil(memory / 13)]) # always require at least two CPUs

  command <<<
    set -euo pipefail
    markers_lines=~{write_lines(n_markers_array)}
    l0_files=~{write_lines(l0_files)}
    echo "~{sep("\n", select_all([bed_file, bim_file, fam_file, covar_file, pheno_file]))}" | \
      cat - $l0_files | tr '\n' '\0' | xargs -0 mv -t .
    cat $l0_files | sed 's/^.*\///;s/.*/& &/;s/[0-9]*$/1/' | while read src dst; do mv --no-clobber $src $dst; done
    sed -i 's/^.*\///;s/[0-9]*$/1/' $l0_files
    awk -F"\t" -v OFS="\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
      {print $(f["FID"]),$(f["IID"]),$(f["~{pheno_name}"])}' "~{basename(pheno_file)}" > "~{basename(pheno_file)}.~{pheno_name}"
    paste $l0_files $markers_lines | \
      awk -F"\t" 'BEGIN {print "~{n_markers} ~{bsize}"}
      {sub("_l0_Y1$","",$1); print $1,1+int(($2-1)/~{bsize}),$2}' > "~{filebase}.master"
    regenie \
      --step 1 \
      --bed "~{basename(bed_file, ".bed")}" \
      --phenoFile "~{basename(pheno_file)}.~{pheno_name}" \
      ~{if defined(covar_file) then "--covarFile \"" + basename(select_first([covar_file])) + "\"" else ""} \
      ~{if binary then "--bt" else ""} \
      --bsize ~{bsize} \
      ~{if loocv then "--loocv" else ""} \
      --run-l1 "~{filebase}.master" \
      --keep-l0 \
      --gz \
      ~{if binary then "--write-null-firth" else ""} \
      ~{if defined(autosome_ct) then "--nauto " + autosome_ct else ""} \
      ~{if defined(regenie_step1_extra_args) then regenie_step1_extra_args else ""} \
      --threads ~{cpu} \
      --out "~{filebase}" \
      1>&2
    rm "~{basename(pheno_file)}.~{pheno_name}" "~{filebase}.master"
    mkdir logs
    mv "~{filebase}.log" "logs/~{filebase}.step1.log"
    mkdir loco
    if [ -f "~{filebase}_1.loco.gz" ]; then
      mv "~{filebase}_1.loco.gz" "loco/~{filebase}.loco.gz"
      sed -i 's/\([^ \t]*\) .*\//\1 /;s/_1\.loco\.gz$/.loco.gz/' "~{filebase}_pred.list"
    else
      touch "loco/~{filebase}.loco.gz"
    fi
    ~{if binary then
      "mkdir firth\n" +
      "if [ -f \"" + filebase + "_1.firth.gz\" ]; then\n" +
      "  mv \"" +  filebase + "_1.firth.gz\" \"firth/" + filebase + ".firth.gz\"\n" +
      "  sed -i 's/\\([^ \\t]*\\) .*\\//\\1 /;s/_1\\.firth\\.gz$/.firth.gz/' \"" + filebase + "_firth.list\"\n" +
      "else\n" +
      "  touch \"firth/" + filebase + ".firth.gz\"\n" +
      "fi" else ""}
    echo "~{sep("\n", select_all([bed_file, bim_file, fam_file, covar_file, pheno_file]))}" | \
      sed 's/^.*\///' | cat - $l0_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File log_file = "logs/" + filebase + ".step1.log"
    String loco_line = read_string(filebase + "_pred.list")
    File loco_file = "loco/" + filebase + ".loco.gz"
    String? firth_line = if binary then read_string(filebase + "_firth.list") else None
    File? firth_file = if binary then "firth/" + filebase + ".firth.gz" else None
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

task serialize_lines {
  input {
    Array[String]+ lines
    String filename

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    lines=~{write_lines(lines)}
    mv $lines "~{filename}"
  >>>

  output {
    File file = filename
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

# number of variants and sample should be input here
# option psam-cols=fid,sex required due to https://github.com/rgcgithub/regenie/issues/105
#
# PLINK2 tries to access more memory than assigned https://groups.google.com/g/plink2-users/c/eLxA5JCwRH0/m/gZmm8RhJAgAJ
# PLINK2 unfortunately needs to read a VCF file twice https://groups.google.com/g/plink2-users/c/hsByNOklyA0/m/ZUHf1MpvAQAJ
# PLINK2 unfortunately cannot split PAR with non-human genomes https://groups.google.com/g/plink2-users/c/88W9O02WXfI
task vcf2pgen {
  input {
    File vcf_file
    File sample_tsv_file
    String dosage_field
    String space_character
    Int? autosome_ct
    Int? par_bp1
    Int? par_bp2

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0

    Float mult = 4.0 # this is to make sure that there is space for the .bcf, the -temporary.pgen, and the .pgen files
    # notice that the .pgen can occupy more space than the .bcf as each encoded dosage requires 2 bytes to be encoded
    # and only dosages that have dosage level identical to their respective genotypes are not encoded
    # see https://docs.juliahub.com/PGENFiles/76R2z/0.1.0/PGEN_description/#Dosages-1
  }

  String filebase = basename(basename(vcf_file, ".bcf"), ".vcf.gz")
  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + mult * vcf_size)])

  command <<<
    set -euo pipefail
    mv "~{vcf_file}" .
    mv "~{sample_tsv_file}" .
    awk 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
      NR>1 {id=$(f["sample_id"]); gsub(" ","~{space_character}",id);
      print 0,id,$(f["computed_gender"])}' "~{basename(sample_tsv_file)}" | \
      sed 's/U$/0/;s/K$/1/' > "~{filebase}.sex"
    bcftools query -l "~{basename(vcf_file)}" | wc -l
    plink2 \
      --threads ~{cpu} \
      --memory ~{round(1024 * memory - 1024)} \
      --bcf "~{basename(vcf_file)}" dosage=~{dosage_field} \
      --update-sex "~{filebase}.sex" \
      --vcf-idspace-to ~{space_character} \
      --const-fid \
      --allow-extra-chr 0 \
      ~{if defined(autosome_ct) then "--chr-set " + autosome_ct else ""} \
      --make-pgen erase-phase psam-cols=fid,sex \
      ~{if !defined(autosome_ct) && defined(par_bp1) && defined(par_bp2) then "--split-par " + select_first([par_bp1]) + " " + select_first([par_bp2]) else ""} \
      --out "~{filebase}" \
      1>&2
    grep -v ^# "~{filebase}.pvar" | wc -l > "~{filebase}.nvar"
    rm "~{filebase}.sex"
    rm "~{basename(vcf_file)}"
    rm "~{basename(sample_tsv_file)}"
  >>>

  output {
    Int n_smpls = read_int(stdout())
    Int n_markers = read_int(filebase + ".nvar")
    File pgen_file = filebase + ".pgen"
    File pvar_file = filebase + ".pvar"
    File psam_file = filebase + ".psam"
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

# loco file needs to be pruned due to bug https://github.com/rgcgithub/regenie/issues/199
# write_lines() hack needed due to bug https://github.com/broadinstitute/cromwell/issues/5540
# see https://github.com/MRCIEU/gwas-vcf-specification for VCF output
task regenie_step2 {
  input {
    String chr
    Int n_phenos
    Int n_covars
    Int n_smpls
    Int n_markers
    File fasta_fai
    File pgen_file
    File pvar_file
    File psam_file
    File? covar_file
    File pheno_file
    String? pop
    Array[String]+ regenie_suffix # suffix array passed due to bug https://github.com/broadinstitute/cromwell/issues/5549
    Boolean binary
    File loco_lst
    Array[File]+ loco_files
    File? firth_lst
    Array[File]? firth_files
    Int bsize
    Float? min_info
    Int? autosome_ct
    String? regenie_step2_extra_args

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  String filebase = basename(pgen_file, ".pgen")
  Float pgen_size = size(pgen_file, "GiB")
  Float pvar_size = size(pvar_file, "GiB")
  Float psam_size = size(psam_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + pgen_size + pvar_size + psam_size)])
  # see print_usage_info() in https://github.com/rgcgithub/regenie/blob/master/src/Regenie.cpp
  Float memory = select_first([memory_override, 3.5 + 8.0 * (3.0 * n_phenos + 2.0 * bsize + n_covars + if binary then (3.0 + n_covars) * n_phenos + 0.5 * bsize else 0) * n_smpls / 1024 / 1024 / 1024])
  Int cpu = select_first([cpu_override, 2 * ceil(memory / 13)]) # always require at least two CPUs

  command <<<
    set -euo pipefail
    loco_files=~{write_lines(loco_files)}
    firth_files=~{if defined(firth_files) then write_lines(select_first([firth_files])) else ""}
    echo "~{sep("\n", select_all([fasta_fai, pgen_file, pvar_file, psam_file, covar_file, pheno_file, loco_lst, firth_lst]))}" | \
      cat - $loco_files~{if defined(firth_files) then " $firth_files" else ""} | tr '\n' '\0' | xargs -0 mv -t .
    ~{if defined(firth_lst) then
      "awk 'NR==FNR {x[$1]++} NR>FNR && $1 in x' \"" + basename(select_first([firth_lst])) + "\" \"" + basename(loco_lst) + "\" > \"" + basename(loco_lst) + ".alt\""
      else ""}
    regenie \
      --step 2 \
      --pgen "~{basename(pgen_file, ".pgen")}" \
      --phenoFile "~{basename(pheno_file)}" \
      ~{if defined(covar_file) then "--covarFile \"" + basename(select_first([covar_file])) + "\"" else ""} \
      --pred "~{if defined(firth_lst) then basename(loco_lst) + ".alt" else basename(loco_lst)}" \
      ~{if binary then "--bt" else ""} \
      --bsize ~{bsize} \
      ~{if binary then "--firth --approx --firth-se" else ""} \
      ~{if binary && defined(firth_lst) then "--use-null-firth \"" +  basename(select_first([firth_lst])) + "\"" else ""} \
      ~{if defined(min_info) then "--minINFO " + min_info else ""} \
      ~{if binary then "--af-cc" else ""} \
      ~{if defined(autosome_ct) then "--nauto " + autosome_ct else ""} \
      ~{if defined(regenie_step2_extra_args) then regenie_step2_extra_args else ""} \
      --gz \
      --threads ~{cpu} \
      --out "~{filebase}" \
      1>&2
    mkdir logs
    mv "~{filebase}.log" "logs/~{filebase}.step2.log"
    ~{if defined(firth_lst) then "rm \"" + basename(loco_lst) + ".alt\"" else ""}

    cut -d" " -f1 "~{basename(loco_lst)}" | while IFS=$'\t' read pheno; do
      bcftools +munge \
        --no-version \
        --columns REGENIE \
        --fai "~{basename(fasta_fai)}" \
        --sample-name $pheno~{if defined(pop) then "." + select_first([pop]) else ""} \
        --output-type b \
        --output "~{filebase + if length(regenie_suffix) > 1 then ".$pheno" else ""}~{if defined(pop) then "." + select_first([pop]) else ""}.gwas.bcf" \
        ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
        "~{filebase}_$pheno.regenie.gz"
      ~{if defined(pop) then "mv \"" + filebase + "_$pheno.regenie.gz\" \"" + filebase + "_$pheno." + select_first([pop]) + ".regenie.gz\"" else ""}
    done
    ~{if length(regenie_suffix) > 1 then
      "awk '{print \"" +  filebase + ".\"$1\"" + (if defined(pop) then "." + select_first([pop]) else "") + ".gwas.bcf\"}' \"" + basename(loco_lst) + "\" | \\\n" +
      "bcftools merge \\\n" +
      "  --no-version \\\n" +
      "  --file-list /dev/stdin \\\n" +
      "  --merge none \\\n" +
      "  --no-index \\\n" +
      "  --output-type b \\\n" +
      "  --output \"" + filebase + (if defined(pop) then "." + select_first([pop]) else "") + ".gwas.bcf\"" +
      (if cpu > 1 then " \\\n  --threads " + (cpu - 1) + "\n" else "\n") +
      "awk '{print \"" +  filebase + ".\"$1\"" + (if defined(pop) then "." + select_first([pop]) else "") + ".gwas.bcf\"}' \"" + basename(loco_lst) + "\" | tr '\\n' '\\0' | xargs -0 rm"
    else ""}
    echo "~{sep("\n", select_all([fasta_fai, pgen_file, pvar_file, psam_file, covar_file, pheno_file, loco_lst, firth_lst]))}" | \
      cat - $loco_files~{if defined(firth_files) then " $firth_files" else ""} | sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File vcf_file = filebase + (if defined(pop) then "." + select_first([pop]) else "") + ".gwas.bcf"
    File log_file = "logs/" + filebase + ".step2.log"
    Array[File] regenie_files = prefix(filebase + "_", regenie_suffix) # suffix array passed due to bug https://github.com/broadinstitute/cromwell/issues/5549
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

# as the VCF to be concatenated can be pretty large and this is only run once, we avoid preemptible computing
task vcf_concat {
  input {
    Array[File]+ vcf_files
    File? ref_fasta
    File? fasta_fai
    File? gff3_file
    File? rsid_vcf_file
    File? rsid_vcf_idx
    String filebase

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 0
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_files, "GiB")
  Float ref_size = size(ref_fasta, "GiB")
  Float dbsnp_size = size(rsid_vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size + ref_size + dbsnp_size)])

  command <<<
    set -euo pipefail
    ~{if defined(ref_fasta) then "mv \"" + select_first([ref_fasta]) + "\" ." else ""}
    ~{if defined(fasta_fai) then "mv \"" + select_first([fasta_fai]) + "\" ." else ""}
    ~{if defined(gff3_file) then "mv \"" + select_first([gff3_file]) + "\" ." else ""}
    ~{if defined(rsid_vcf_file) then "mv \"" + select_first([rsid_vcf_file]) + "\" ." else ""}
    ~{if defined(rsid_vcf_idx) then "mv \"" + select_first([rsid_vcf_idx]) + "\" ." else ""}
    vcf_files=~{write_lines(vcf_files)}
    cat $vcf_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $vcf_files
    bcftools concat \
      --no-version \
      --output-type ~{if defined(gff3_file) then "u" else "b"} \
      --file-list $vcf_files \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
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
    ~{if defined(ref_fasta) then "rm \"" + basename(select_first([ref_fasta])) + "\"" else ""}
    ~{if defined(fasta_fai) then "rm \"" + basename(select_first([fasta_fai])) + "\"" else ""}
    ~{if defined(gff3_file) then "rm \"" + basename(select_first([gff3_file])) + "\"" else ""}
    cat $vcf_files | tr '\n' '\0' | xargs -0 rm
    ~{if defined(rsid_vcf_file) then
      "bcftools annotate \\\n" +
      "  --no-version \\\n" +
      "  --annotations \"" + basename(select_first([rsid_vcf_file])) + "\" \\\n" +
      "  --columns RS \\\n" +
      "  --output \"" + filebase + ".bcf\" \\\n" +
      "  --output-type b \\\n" +
      "  \"" + filebase + ".tmp.bcf\" \\\n" +
      "  --write-index\n" +
      "  rm \"" + filebase + ".tmp.bcf\" \"" + filebase + ".tmp.bcf.csi\" \"" + basename(select_first([rsid_vcf_file])) + "\""
      else ""}
    ~{if defined(rsid_vcf_idx) then "rm \"" + basename(select_first([rsid_vcf_idx])) + "\"" else ""}
  >>>

  output {
    File vcf_file = filebase + ".bcf"
    File vcf_idx = filebase + ".bcf.csi"
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

task assoc_concat {
  input {
    Array[File]+ assoc_files
    Int n_x_chr
    Boolean zst = false
    String filebase

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float assoc_size = size(assoc_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * assoc_size)])

  command <<<
    set -euo pipefail
    assoc_files=~{write_lines(assoc_files)}
    cat $assoc_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $assoc_files
    cat $assoc_files | tr '\n' '\0' | xargs -0 ~{if zst then "-n1 plink2 --zst-decompress" else "zcat"} | \
      awk 'NR==1 || $0!~"^CHROM" && $0!~"^#CHROM"' | \
      sed 's/^CHROM/#CHROM/;s/^~{n_x_chr}/X/;s/^PAR[12]/X/' | tr ' ' '\t' | \
      bgzip > "~{filebase}"
    tabix --begin 2 --end 2 --force "~{filebase}"
    tabix --only-header "~{filebase}" | wc -l
    cat $assoc_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File file = filebase
    File index = filebase + ".tbi"
    Boolean has_data = if read_int(stdout()) == 0 then false else true
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

task assoc_plot {
  input {
    File assoc_file
    String? genome
    File? cyto_file
    Int? autosome_ct
    Float? min_af
    String filebase

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float assoc_size = size(assoc_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + assoc_size)])

  command <<<
    set -euo pipefail
    mv "~{assoc_file}" .
    ~{if defined(cyto_file) then "mv \"" + select_first([cyto_file]) + "\" ." else ""}
    assoc_plot.R \
      ~{if defined(cyto_file) then "--cytoband \"" + basename(select_first([cyto_file])) + "\"" else
        if defined(genome) then "--genome \"" + select_first([genome]) + "\"" else ""} \
      ~{if defined(autosome_ct) then "--nauto " + select_first([autosome_ct]) else ""} \
      ~{if defined(min_af) then "--min-af " + select_first([min_af]) else ""} \
      --tbx "~{basename(assoc_file)}" \
      --png "~{filebase}.png"
    rm "~{basename(assoc_file)}"
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

# use of !(a!=b) due to bug Cromwell team will not fix: https://github.com/broadinstitute/cromwell/issues/5602
# plink2 will return error code 7 when covariate-only Firth regression fails to converge or when variance inflation factor is too high
# plink2 will return error code 13 when no diploid variants remain for --glm hetonly
# when PLINK2 updates to version Mar 18 2024 or newer, you can remove the single-prec-cc otions https://groups.google.com/g/plink2-users/c/4oKzPt_Xi74/m/z--hBJV-AgAJ
task plink_glm {
  input {
    String chr_num
    String pheno_name
    Int n_phenos
    Int n_covars
    Int n_smpls
    Int n_markers
    File pgen_file
    File pvar_file
    File psam_file
    File? loco_file
    File? covar_file
    File pheno_file
    Boolean binary
    Int? max_vif
    Float? max_corr
    String? plink_extra_args

    String docker
    Int cpu = 1
    Int? disk_size_override
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  String filebase = basename(pgen_file, ".pgen")
  Float pgen_size = size(pgen_file, "GiB")
  Float pvar_size = size(pvar_file, "GiB")
  Float psam_size = size(psam_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + pgen_size + pvar_size + psam_size)])

  command <<<
    set -euo pipefail
    echo "~{sep("\n", select_all([pgen_file, pvar_file, psam_file, covar_file, pheno_file, loco_file]))}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    ~{if defined(loco_file) then
      "zcat \"" + basename(select_first([loco_file])) + "\" | \\\n" +
      "awk 'BEGIN {print \"IID\\tLOCO\"} NR==1 {for (i=2; i<=NF; i++) {sub(\"^0_\", \"\", $i); f[i] = $i}}\n" +
      "  $1==\"" + chr_num + "\" {for (i=2; i<=NF; i++) print f[i]\"\\t\"$i}' " +
      (if defined(covar_file) then "| \\\n" +
        "awk -F\"\\t\" -v OFS=\"\\t\" 'NR==FNR {x[$1]=$2} NR>FNR && $2 in x {print $0\"\\t\"x[$2]}' - \"" + basename(select_first([covar_file])) + "\" "
        else "") +
      "> \"" +  filebase + ".cov\""
    else ""}
    plink2 \
      --threads ~{cpu} \
      --memory ~{round(1024 * memory - 512)} \
      --pgen "~{basename(pgen_file)}" \
      --pvar "~{basename(pvar_file)}" \
      --psam "~{basename(psam_file)}" \
      ~{if defined(loco_file) then "--covar \"" + filebase + ".cov\"" else
        if defined(covar_file) then "--covar \"" + basename(select_first([covar_file])) + "\"" else ""} \
      ~{if defined(loco_file) || defined(covar_file) then "--covar-variance-standardize" else ""} \
      ~{if defined(loco_file) || defined(covar_file) then "--require-covar" else ""} \
      --1 --pheno "~{basename(pheno_file)}" \
      --pheno-name ~{pheno_name} \
      --require-pheno \
      --glm zs log10 hetonly hide-covar~{if !defined(loco_file) && !defined(covar_file) then " allow-no-covars" else ""} cc-residualize single-prec-cc cols=+a1freq,+machr2 \
      ~{if defined(max_vif) then "--vif " + select_first([max_vif]) else ""} \
      ~{if defined(max_corr) then "--max-corr " + select_first([max_corr]) else ""} \
      ~{if defined(plink_extra_args) then plink_extra_args else ""} \
      --out "~{filebase}" \
      1>&2 || if [[ $? -eq 7 || $? -eq 13 ]]; then
        echo -en "\x28\xb5\x2f\xfd\x24\x00\x01\x00\x00\x99\xe9\xd8\x51" > "~{filebase}.~{pheno_name}.glm.~{if binary then "logistic.hybrid" else "linear"}.zst"
      else exit $?; fi
    ~{if defined(loco_file) then "rm \"" + filebase + ".cov\"" else ""}
    mkdir logs
    mv "~{filebase}.log" "logs/~{filebase}.~{pheno_name}.glm.log"
    echo "~{sep("\n", select_all([pgen_file, pvar_file, psam_file, covar_file, pheno_file, loco_file]))}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File log_file = "logs/" + filebase + "." + pheno_name + ".glm.log"
    File assoc_file = filebase + "." + pheno_name + ".glm." + (if binary then "logistic.hybrid" else "linear") + ".zst"
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
