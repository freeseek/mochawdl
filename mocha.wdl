version development

## Copyright (c) 2020-2024 Giulio Genovese
##
## Version 2024-05-05
##
## Contact Giulio Genovese <giulio.genovese@gmail.com>
##
## This WDL workflow runs MoChA on a cohort of samples genotyped with either Illumina or Affymetrix DNA microarrays
##
## Cromwell version support
## - Successfully tested on v86
##
## Distributed under terms of the MIT License

struct Reference {
  String name
  File fasta
  Int min_chr_len
  File? rules_file
  String? mhc_reg
  String? kir_reg
  String? nonpar_reg
  File? dup_file
  File genetic_map_file
  File? cnp_file
  File cyto_file
  String? panel_pfx
  String? panel_sfx
  String panel_idx
  Int? n_panel_smpls
}

workflow mocha {
  input {
    String sample_set_id
    String mode # idat gtc cel chp txt vcf pvcf
    String target = "calls" # vcf pvcf calls pngs
    Boolean realign = false
    Boolean wgs = false
    Boolean gtc_output = false # only for idat mode
    Boolean chp_output = false # only for cel mode
    Int idat_batch_size = 256
    Int gtc_batch_size = 1024
    Int chp_batch_size = 1024
    Float max_win_size_cm = 50.0
    Float overlap_size_cm = 5.0
    Float sample_call_rate_thr = 0.97
    Float variant_call_rate_thr = 0.97
    Float baf_auto_thr = 0.03
    String ext_string = "as"

    String ref_name = "GRCh38"
    String? ref_path
    String? ref_fasta
    Int? min_chr_len
    File? ref_rules_file
    String? mhc_reg
    String? kir_reg
    String? nonpar_reg
    String? dup_file
    String? genetic_map_file
    String? cnp_file
    String? cyto_file
    String? panel_pfx
    String? panel_sfx
    String? panel_idx
    Int? n_panel_smpls

    String manifest_path = ""
    File sample_tsv_file # sample_id batch_id green_idat red_idat gtc cel chp computed_gender call_rate
    File batch_tsv_file # batch_id path bpm csv egt sam xml zip probeset_ids snp report calls confidences summary vcf vcf_index xcl_vcf xcl_vcf_index
    String? data_path
    File? pedigree_file
    File? duplicate_samples_file
    File? extra_xcl_vcf_file
    String? gtc2vcf_extra_args
    String? phase_extra_args
    String? mocha_extra_args
    String? mocha_plot_extra_args
    String basic_bash_docker = "debian:stable-slim"
    String pandas_docker = "amancevice/pandas:slim"
    String docker_repository = "us.gcr.io/mccarroll-mocha"
    String bcftools_docker = "bcftools:1.20-20240505"
    String apt_docker = "apt:1.20-20240505"
    String shapeit5_docker = "shapeit5:1.20-20240505"
    String r_mocha_docker = "r_mocha:1.20-20240505"
    Boolean? table_output
    Boolean do_not_use_reference = false
    Boolean use_shapeit5_ligate = false
    String delim = "~"
    Array[String]? chip_type
    String tags = "GT,BAF,LRR"
    Int? gc_window_size
  }

  Boolean mode_is_vcf = mode == "vcf" || mode == "pvcf" || wgs

  String docker_repository_with_sep = docker_repository + if docker_repository != "" && docker_repository == sub(docker_repository, "/$", "") then "/" else ""

  String ref_path_with_sep = select_first([ref_path, ""]) + if defined(ref_path) && select_first([ref_path]) == sub(select_first([ref_path]), "/$", "") then "/" else ""
  Reference ref = object {
    name: ref_name,
    fasta: if defined(ref_fasta) then ref_path_with_sep + select_first([ref_fasta]) else if ref_name == "GRCh38" then ref_path_with_sep + "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" else if ref_name == "GRCh37" then ref_path_with_sep + "human_g1k_v37.fasta" else None,
    min_chr_len: select_first([min_chr_len, 2000000]),
    rules_file: if defined(ref_rules_file) then ref_path_with_sep + select_first([ref_rules_file]) else None,
    mhc_reg: if defined(mhc_reg) then select_first([mhc_reg]) else if ref_name == "GRCh38" then "chr6:27518932-33480487" else if ref_name == "GRCh37" then "6:27486711-33448264" else None,
    kir_reg: if defined(kir_reg) then select_first([kir_reg]) else if ref_name == "GRCh38" then "chr19:54071493-54992731" else if ref_name == "GRCh37" then "19:54574747-55504099" else None,
    nonpar_reg: if defined(nonpar_reg) then select_first([nonpar_reg]) else if ref_name == "GRCh38" then "chrX:2781479-155700628" else if ref_name == "GRCh37" then "X:2699520-154930289" else None,
    dup_file: if defined(dup_file) then ref_path_with_sep + select_first([dup_file]) else if ref_name == "GRCh38" || ref_name == "GRCh37" then ref_path_with_sep + "segdups.bed.gz" else None,
    genetic_map_file: if defined(genetic_map_file) then ref_path_with_sep + select_first([genetic_map_file]) else if ref_name == "GRCh38" then ref_path_with_sep + "genetic_map_hg38_withX.txt.gz" else if ref_name == "GRCh37" then ref_path_with_sep + "genetic_map_hg19_withX.txt.gz" else None,
    cnp_file: if defined(cnp_file) then ref_path_with_sep + select_first([cnp_file]) else if ref_name == "GRCh38" || ref_name == "GRCh37" then ref_path_with_sep + "cnps.bed" else None,
    cyto_file: if defined(cyto_file) then ref_path_with_sep + select_first([cyto_file]) else if ref_name == "GRCh38" || ref_name == "GRCh37" then ref_path_with_sep + "cytoBand.txt.gz" else None,
    panel_pfx: if defined(panel_pfx) then ref_path_with_sep + select_first([panel_pfx]) else if ref_name == "GRCh38" then ref_path_with_sep + "1kGP_high_coverage_Illumina." else if ref_name == "GRCh37" then ref_path_with_sep + "ALL.chr" else None,
    panel_sfx: if defined(panel_sfx) then select_first([panel_sfx]) else if ref_name == "GRCh38" then ".bcf" else if ref_name == "GRCh37" then ".phase3_integrated.20130502.genotypes.bcf" else None,
    panel_idx: select_first([panel_idx, ".csi"]),
    n_panel_smpls: if defined(n_panel_smpls) then select_first([n_panel_smpls]) else if ref_name == "GRCh38" then 3202 else if ref_name == "GRCh37" then 2504 else None,
  }

  # read table with batches information (scatter could be avoided if there was a tail() function)
  call tsv_sorted as batch_sorted_tsv { input: tsv_file = batch_tsv_file, column = "batch_id", check_dups = true, docker = basic_bash_docker }
  Array[Array[String]] batch_tsv = read_tsv(batch_sorted_tsv.file)
  Int n_batches = length(batch_tsv)-1
  scatter (idx in range(n_batches)) { Array[String] batch_tsv_rows = batch_tsv[(idx+1)] }
  Map[String, Array[String]] batch_tbl = as_map(zip(batch_tsv[0], transpose(batch_tsv_rows)))
  Array[String] batches = batch_tbl["batch_id"]
  # check which variables are in batch table (see https://github.com/openwdl/wdl/issues/305)
  Boolean is_path_in_batch_tbl = length(collect_by_key(zip(flatten([keys(batch_tbl),["path"]]),range(length(keys(batch_tbl))+1)))["path"])>1
  Boolean is_sam_in_batch_tbl = length(collect_by_key(zip(flatten([keys(batch_tbl),["sam"]]),range(length(keys(batch_tbl))+1)))["sam"])>1
  Boolean is_probeset_ids_in_batch_tbl = length(collect_by_key(zip(flatten([keys(batch_tbl),["probeset_ids"]]),range(length(keys(batch_tbl))+1)))["probeset_ids"])>1
  Boolean is_confidences_in_batch_tbl = length(collect_by_key(zip(flatten([keys(batch_tbl),["confidences"]]),range(length(keys(batch_tbl))+1)))["confidences"])>1

  String manifest_path_with_sep = manifest_path + (if manifest_path == "" || sub(manifest_path, "/$", "") != manifest_path then "" else "/")
  scatter (idx in range(n_batches)) {
    String data_paths_with_sep = (if defined(data_path) then sub(select_first([data_path]), "/$", "") + "/" else "") +
                                 (if is_path_in_batch_tbl then sub(batch_tbl["path"][idx], "/$", "") + "/" else "")
    String pfxs = sample_set_id + (if n_batches == 1 then "" else "." + batches[idx])
  }

  # aligns manifest file to human genome reference if requested
  if (realign && !mode_is_vcf) {
    # hack due to lack of unique function
    Array[String] csv_files = keys(collect_by_key(zip(batch_tbl["csv"], batches)))
    scatter (csv_file in csv_files) {
      call csv2bam {
        input:
          plugin = if mode == "idat" || mode == "gtc" then "gtc2vcf" else "affy2vcf",
          csv_file = manifest_path_with_sep + csv_file,
          ref_fasta = ref.fasta,
          ref_fasta_idxs = prefix(ref.fasta + ".", ["amb", "ann", "bwt", "pac", "sa"]),
          docker = docker_repository_with_sep + bcftools_docker,
      }
    }
    Map[String, File] csv2sam = as_map(zip(csv_files, csv2bam.bam_file))
  }
  scatter (idx in range(n_batches)) {
    File? sams = if realign && !mode_is_vcf then select_first([csv2sam])[(batch_tbl["csv"][idx])]
      else if is_sam_in_batch_tbl then manifest_path_with_sep + batch_tbl["sam"][idx] else None
  }

  Array[Array[String]] ref_fasta_fai_tbl = transpose(read_tsv(ref.fasta + ".fai"))
  scatter (idx in range(length(ref_fasta_fai_tbl[0]))) {
    Int fai_len = ref_fasta_fai_tbl[1][idx]
    if (fai_len > ref.min_chr_len && ref_fasta_fai_tbl[0][idx] != "Y" && ref_fasta_fai_tbl[0][idx] != "chrY") {
      String chrs = ref_fasta_fai_tbl[0][idx]
      Int lens = fai_len
    }
  }

  if (!mode_is_vcf) {
    # resort table with sample information and extract sample_id column
    call tsv_sorted as sample_sorted_tsv { input: tsv_file = sample_tsv_file, column = "batch_id", docker = basic_bash_docker }
    call tsv_column as sample_id_lines { input: tsv_file = sample_sorted_tsv.file, column = "sample_id", check_dups = true, docker = basic_bash_docker }
    call tsv_column as batch_id_lines { input: tsv_file = sample_sorted_tsv.file, column = "batch_id", docker = basic_bash_docker }

    # process Illumina data
    if (mode == "idat") {
      call tsv_column as green_idat_lines { input: tsv_file = sample_sorted_tsv.file, column = "green_idat", docker = basic_bash_docker }
      call tsv_column as red_idat_lines { input: tsv_file = sample_sorted_tsv.file, column = "red_idat", docker = basic_bash_docker }
      # group samples by IDAT batches
      call batch_scatter as idat { input: batch_id_file = batch_id_lines.file, sub_batch_size = idat_batch_size, delim = delim, docker = basic_bash_docker }
      Map[String, Array[String]] idat_batch2green_idat_files = collect_by_key(zip(read_lines(idat.sub_batch_id), read_lines(green_idat_lines.file)))
      Map[String, Array[String]] idat_batch2red_idat_files = collect_by_key(zip(read_lines(idat.sub_batch_id), read_lines(red_idat_lines.file)))
      scatter (idx in range(length(idat.sub_batches))) {
        call idat2gtc {
          input:
            bpm_file = manifest_path_with_sep + batch_tbl["bpm"][(idat.idxs[idx])],
            egt_file = manifest_path_with_sep + batch_tbl["egt"][(idat.idxs[idx])],
            green_idat_files = prefix(data_paths_with_sep[(idat.idxs[idx])], idat_batch2green_idat_files[(idat.sub_batches[idx])]),
            red_idat_files = prefix(data_paths_with_sep[(idat.idxs[idx])], idat_batch2red_idat_files[(idat.sub_batches[idx])]),
            filebase = sample_set_id + "." + idat.sub_batches[idx],
            docker = docker_repository_with_sep + bcftools_docker,
        }
      }
      call tsv_concat as green_idat_tsv { input: tsv_files = idat2gtc.green_idat_tsv, filebase = sample_set_id + ".green_idat", docker = basic_bash_docker }
      call tsv_concat as red_idat_tsv { input: tsv_files = idat2gtc.red_idat_tsv, filebase = sample_set_id + ".red_idat", docker = basic_bash_docker }
    }

    if (mode == "gtc") {
      call tsv_column as gtc_lines { input: tsv_file = sample_sorted_tsv.file, column = "gtc", docker = basic_bash_docker }
    }

    if (mode == "idat" || mode == "gtc") {
      # group samples by GTC batches
      call batch_scatter as gtc { input: batch_id_file = batch_id_lines.file, sub_batch_size = gtc_batch_size, delim = delim, docker = basic_bash_docker }
      call get_reheader_maps as gtc_reheader { input: batches = gtc.sub_batches, barcodes_file = select_first([green_idat_lines.file, gtc_lines.file]), sample_id_file = sample_id_lines.file, batch_id_file = gtc.sub_batch_id, docker = basic_bash_docker }
      Array[String]+ input_gtc_files = if mode == "idat" then flatten(select_first([idat2gtc.gtc_files])) else read_lines(select_first([gtc_lines.file]))
      Map[String, Array[String]] gtc_batch2gtc_files = collect_by_key(zip(read_lines(gtc.sub_batch_id), input_gtc_files))
      scatter (idx in range(length(gtc.sub_batches))) {
        call gtc2vcf {
          input:
            tags = tags,
            bpm_file = manifest_path_with_sep + batch_tbl["bpm"][(gtc.idxs[idx])],
            csv_file = manifest_path_with_sep + batch_tbl["csv"][(gtc.idxs[idx])],
            egt_file = manifest_path_with_sep + batch_tbl["egt"][(gtc.idxs[idx])],
            ref_fasta = ref.fasta,
            ref_fasta_fai = ref.fasta + ".fai",
            gc_window_size = gc_window_size,
            gtc_files = prefix(if mode == "idat" then "" else data_paths_with_sep[(gtc.idxs[idx])], gtc_batch2gtc_files[(gtc.sub_batches[idx])]),
            gtc2vcf_extra_args = gtc2vcf_extra_args,
            sam_file = sams[(gtc.idxs[idx])],
            reheader_file = gtc_reheader.map_files[idx],
            filebase = sample_set_id + "." + gtc.sub_batches[idx],
            docker = docker_repository_with_sep + bcftools_docker,
        }
      }
      call tsv_concat as gtc_tsv { input: tsv_files = gtc2vcf.gtc_tsv, filebase = sample_set_id + ".gtc", docker = basic_bash_docker }

      # this job can be long, so it is better to run as non-preemptible
      Map[Int, Array[String]] idx2gtc2vcf_files = collect_by_key(zip(gtc.idxs, gtc2vcf.vcf_file))
      Map[Int, Array[String]] idx2gtc2vcf_idxs = collect_by_key(zip(gtc.idxs, gtc2vcf.vcf_idx))
      Map[Int, Array[String]] idx2gtc2vcf_n_smpls = collect_by_key(zip(gtc.idxs, gtc2vcf.n_smpls))
      scatter (idx in range(n_batches)) {
        if (length(idx2gtc2vcf_files[idx]) > 1) {
          call vcf_merge as gtc2vcf_merge {
            input:
              vcf_files = idx2gtc2vcf_files[idx],
              filebase = pfxs[idx],
              docker = docker_repository_with_sep + bcftools_docker
          }
        }
        File gtc2vcf_files = select_first([gtc2vcf_merge.vcf_file, idx2gtc2vcf_files[idx][0]])
        File gtc2vcf_idxs = select_first([gtc2vcf_merge.vcf_idx, idx2gtc2vcf_idxs[idx][0]])
        Int gtc2vcf_n_smpls = select_first([gtc2vcf_merge.n_smpls, idx2gtc2vcf_n_smpls[idx][0]])
      }
    }

    if (mode == "cel" || mode == "txt") {
      call tsv_column as cel_lines { input: tsv_file = sample_sorted_tsv.file, column = "cel", docker = basic_bash_docker }
    }

    # process Affymetrix data
    if (mode == "cel") {
      Map[String, Array[String]] batch2cel_files = collect_by_key(zip(read_lines(batch_id_lines.file), read_lines(select_first([cel_lines.file]))))
      scatter (idx in range(n_batches)) {
        call cel2affy as cel2chp {
          input:
            xml_file = manifest_path_with_sep + batch_tbl["xml"][idx],
            zip_file = manifest_path_with_sep + batch_tbl["zip"][idx],
            cel_files = prefix(data_paths_with_sep[idx], batch2cel_files[(batches[idx])]),
            probeset_file = if is_probeset_ids_in_batch_tbl && batch_tbl["probeset_ids"][idx] != "" then manifest_path_with_sep + batch_tbl["probeset_ids"][idx] else None,
            chip_type = chip_type,
            table_output = table_output,
            filebase = pfxs[idx],
            docker = docker_repository_with_sep + apt_docker
        }
      }
      call tsv_concat as cel_tsv { input: tsv_files = cel2chp.cel_tsv, filebase = sample_set_id + ".cel", docker = basic_bash_docker }
    }

    if (mode == "chp") {
      call tsv_column as chp_lines { input: tsv_file = sample_sorted_tsv.file, column = "chp", docker = basic_bash_docker }
    }

    if (mode == "cel" || mode == "chp") {
      # group samples by CHP batches
      call batch_scatter as chp { input: batch_id_file = batch_id_lines.file, sub_batch_size = chp_batch_size, delim = delim, docker = basic_bash_docker }
      call get_reheader_maps as chp_reheader { input: batches = chp.sub_batches, barcodes_file = select_first([cel_lines.file, chp_lines.file]), sample_id_file = sample_id_lines.file, batch_id_file = chp.sub_batch_id, docker = basic_bash_docker }
      Array[String]+ input_chp_files = if mode == "cel" then flatten(select_first([cel2chp.chp_files])) else read_lines(select_first([chp_lines.file]))
      Map[String, Array[String]] chp_batch2chp_files = collect_by_key(zip(read_lines(chp.sub_batch_id), input_chp_files))
      scatter (idx in range(length(chp.sub_batches))) {
        call chp2vcf {
          input:
            tags = tags,
            csv_file = manifest_path_with_sep + batch_tbl["csv"][(chp.idxs[idx])],
            ref_fasta = ref.fasta,
            ref_fasta_fai = ref.fasta + ".fai",
            gc_window_size = gc_window_size,
            probeset_file = if is_probeset_ids_in_batch_tbl && batch_tbl["probeset_ids"][(chp.idxs[idx])] != "" then manifest_path_with_sep + batch_tbl["probeset_ids"][(chp.idxs[idx])] else None,
            snp_file = if mode == "cel" then select_first([cel2chp.snp_file])[(chp.idxs[idx])] else data_paths_with_sep[(chp.idxs[idx])] + batch_tbl["snp"][(chp.idxs[idx])],
            chp_files = prefix(if mode == "cel" then "" else data_paths_with_sep[(chp.idxs[idx])], chp_batch2chp_files[(chp.sub_batches[idx])]),
            sam_file = sams[(chp.idxs[idx])],
            reheader_file = chp_reheader.map_files[idx],
            filebase = sample_set_id + "." + chp.sub_batches[idx],
            docker = docker_repository_with_sep + bcftools_docker
        }
      }

      # this job can be long, so it is better to run as non-preemptible
      Map[Int, Array[String]] idx2chp2vcf_files = collect_by_key(zip(chp.idxs, chp2vcf.vcf_file))
      Map[Int, Array[String]] idx2chp2vcf_idxs = collect_by_key(zip(chp.idxs, chp2vcf.vcf_idx))
      Map[Int, Array[String]] idx2chp2vcf_n_smpls = collect_by_key(zip(chp.idxs, chp2vcf.n_smpls))
      scatter (idx in range(n_batches)) {
        if (length(idx2chp2vcf_files[idx]) > 1) {
          call vcf_merge as chp2vcf_merge {
            input:
              vcf_files = idx2chp2vcf_files[idx],
              filebase = pfxs[idx],
              docker = docker_repository_with_sep + bcftools_docker
          }
        }
        File chp2vcf_files = select_first([chp2vcf_merge.vcf_file, idx2chp2vcf_files[idx][0]])
        File chp2vcf_idxs = select_first([chp2vcf_merge.vcf_idx, idx2chp2vcf_idxs[idx][0]])
        Int chp2vcf_n_smpls = select_first([chp2vcf_merge.n_smpls, idx2chp2vcf_n_smpls[idx][0]])
      }
    }

    if (mode == "txt") {
      call get_reheader_maps as txt_reheader { input: batches = batches, barcodes_file = select_first([cel_lines.file]), sample_id_file = sample_id_lines.file, batch_id_file = batch_id_lines.file, docker = basic_bash_docker }
      scatter (idx in range(n_batches)) {
        # this job can be long, so it is better to run as non-preemptible
        call txt2vcf {
          input:
            tags = tags,
            csv_file = manifest_path_with_sep + batch_tbl["csv"][idx],
            ref_fasta = ref.fasta,
            ref_fasta_fai = ref.fasta + ".fai",
            gc_window_size = gc_window_size,
            probeset_file = if is_probeset_ids_in_batch_tbl && batch_tbl["probeset_ids"][idx] != ""
              then manifest_path_with_sep + batch_tbl["probeset_ids"][idx] else None,
            calls_file = data_paths_with_sep[idx] + batch_tbl["calls"][idx],
            confidences_file = if is_confidences_in_batch_tbl then data_paths_with_sep[idx] + batch_tbl["confidences"][idx] else None,
            summary_file = data_paths_with_sep[idx] + batch_tbl["summary"][idx],
            report_file = data_paths_with_sep[idx] + batch_tbl["report"][idx],
            snp_file = data_paths_with_sep[idx] + batch_tbl["snp"][idx],
            sam_file = sams[idx],
            reheader_file = txt_reheader.map_files[idx],
            filebase = pfxs[idx],
            docker = docker_repository_with_sep + bcftools_docker
        }
      }
    }

    if (mode == "cel" || mode == "chp" || mode == "txt") {
      call tsv_concat as affy_tsv { input: tsv_files = select_first([chp2vcf.affy_tsv, txt2vcf.affy_tsv]), filebase = sample_set_id + ".affy", docker = basic_bash_docker }
    }

    Array[File] unphased_vcf_files = select_first([gtc2vcf_files, chp2vcf_files, txt2vcf.vcf_file])
    Array[File] unphased_vcf_idxs = select_first([gtc2vcf_idxs, chp2vcf_idxs, txt2vcf.vcf_idx])
    Array[Int] unphased_vcf_n_smpls = select_first([gtc2vcf_n_smpls, chp2vcf_n_smpls, txt2vcf.n_smpls])
    call lst_flatten as flatten_sample_id_lines { input: lst_files = select_first([gtc2vcf.sample_id_lines, chp2vcf.sample_id_lines, txt2vcf.sample_id_lines]), filebase = "sample_id", docker = basic_bash_docker }
    call tsv_column as override_computed_gender_lines { input: tsv_file = sample_sorted_tsv.file, column = "computed_gender", fail = false, docker = basic_bash_docker }
    if (override_computed_gender_lines.failed) {
      call tsv_column as computed_gender_lines { input: tsv_file = select_first([gtc_tsv.file, affy_tsv.file, sample_tsv_file]), column = "computed_gender", docker = basic_bash_docker }
    }
    call tsv_column as override_call_rate_lines { input: tsv_file = sample_sorted_tsv.file, column = "call_rate", fail = false, docker = basic_bash_docker }
    if (override_call_rate_lines.failed) {
      call tsv_column as call_rate_lines { input: tsv_file = select_first([gtc_tsv.file, affy_tsv.file, sample_tsv_file]), column = "call_rate", docker = basic_bash_docker }
    }
    call lst_paste as sample_tsv {
      input:
        lst_files = select_all([flatten_sample_id_lines.file, select_first([computed_gender_lines.file, override_computed_gender_lines.file]), select_first([call_rate_lines.file, override_call_rate_lines.file])]),
        headers = ['sample_id', 'computed_gender', 'call_rate'],
        filebase = sample_set_id + ".sample",
        docker = basic_bash_docker
    }
  }

  if (mode != "pvcf" && target != "vcf") {
    call ref_scatter {
      input:
        chrs = select_all(chrs),
        lens = select_all(lens),
        genetic_map_file = ref.genetic_map_file,
        max_win_size_cm = max_win_size_cm,
        overlap_size_cm = overlap_size_cm,
        docker = pandas_docker
    }
    scatter (idx in range(n_batches)) {
      call vcf_scatter {
        input:
          vcf_file = if mode_is_vcf then data_paths_with_sep[idx] + batch_tbl["vcf"][idx] else select_first([unphased_vcf_files])[idx],
          intervals_tsv = ref_scatter.intervals_tsv,
          docker = docker_repository_with_sep + bcftools_docker
      }
    }
    if (defined(extra_xcl_vcf_file)) {
      call vcf_scatter as xcl_vcf_scatter {
        input:
          vcf_file = select_first([extra_xcl_vcf_file]),
          intervals_tsv = ref_scatter.intervals_tsv,
          docker = docker_repository_with_sep + bcftools_docker
      }
    }

    call lst_concat as sample_id_split_tsv { input: lst_files = vcf_scatter.sample_id_lines, filebase = "split_sample_id", docker = basic_bash_docker }
    Int n_smpls = length(flatten(read_tsv(sample_id_split_tsv.file)))
    Boolean use_reference = if do_not_use_reference then false else n_smpls <= 2 * ref.n_panel_smpls
    Array[Array[File]] interval_slices = transpose(vcf_scatter.vcf_files)
    Array[Array[String]] intervals_tbl = transpose(read_tsv(ref_scatter.intervals_tsv))
    if (defined(ref.mhc_reg)) {
      String? mhc_chr = sub(select_first([ref.mhc_reg]), ":.*$", "")
      Int mhc_beg = sub(sub(select_first([ref.mhc_reg]), "-.*$", ""), "^.*:", "")
      Int mhc_end = sub(select_first([ref.mhc_reg]), "^.*-", "")
    }
    if (defined(ref.nonpar_reg)) {
      String? nonpar_chr = sub(select_first([ref.nonpar_reg]), ":.*$", "")
      Int nonpar_beg = sub(sub(select_first([ref.nonpar_reg]), "-.*$", ""), "^.*:", "")
      Int nonpar_end = sub(select_first([ref.nonpar_reg]), "^.*-", "")
    }
    scatter (idx in range(length(intervals_tbl[0]))) {
      if (length(interval_slices[idx]) > 1) {
        call vcf_merge {
          input:
            vcf_files = interval_slices[idx],
            filebase = sample_set_id + "." + idx,
            docker = docker_repository_with_sep + bcftools_docker
        }
      }
      call vcf_qc {
        input:
          vcf_file = select_first([vcf_merge.vcf_file, interval_slices[idx][0]]),
          vcf_idx = vcf_merge.vcf_idx,
          dup_file = ref.dup_file,
          sample_tsv_file = select_first([sample_tsv.file, sample_tsv_file]),
          sample_call_rate_thr = if wgs then None else sample_call_rate_thr,
          variant_call_rate_thr = variant_call_rate_thr,
          duplicate_samples_file = duplicate_samples_file,
          extra_xcl_vcf_file = if defined(extra_xcl_vcf_file) then select_first([xcl_vcf_scatter.vcf_files])[idx] else None,
          docker = docker_repository_with_sep + bcftools_docker
      }
      Int buffer_beg = intervals_tbl[1][idx] # cast string to integer
      Int buffer_end = intervals_tbl[2][idx] # cast string to integer
      Int beg = intervals_tbl[3][idx] # cast string to integer
      Int end = intervals_tbl[4][idx] # cast string to integer
      call vcf_shapeit5 {
        input:
          n_smpls = n_smpls,
          n_markers = vcf_qc.n_markers,
          unphased_vcf_file = vcf_qc.qc_vcf_file,
          unphased_vcf_idx = vcf_qc.qc_vcf_idx,
          genetic_map_file = ref.genetic_map_file,
          n_chrs = length(select_all(chrs)),
          pedigree_file = pedigree_file,
          sample_tsv_file = if defined(ref.nonpar_reg) && intervals_tbl[0][idx] == select_first([nonpar_chr]) && beg >= select_first([nonpar_beg]) && end <= select_first([nonpar_end]) then select_first([sample_tsv.file, sample_tsv_file]) else None,
          n_panel_smpls = if use_reference then ref.n_panel_smpls else None,
          ref_vcf_file = if use_reference then ref.panel_pfx + intervals_tbl[0][idx] + ref.panel_sfx else None,
          ref_vcf_idx = if use_reference then ref.panel_pfx + intervals_tbl[0][idx] + ref.panel_sfx + ref.panel_idx else None,
          ref_fasta_fai = ref.fasta + ".fai",
          input_region = intervals_tbl[0][idx] + ":" + (1 + beg) + "-" + end,
          scaffold_region = intervals_tbl[0][idx] + ":" + (1 + buffer_beg) + "-" + buffer_end,
          chr = intervals_tbl[0][idx],
          mhc = if defined(ref.mhc_reg) then intervals_tbl[0][idx] == select_first([mhc_chr]) && beg <= select_first([mhc_end]) && end > select_first([mhc_beg]) else false,
          maf = if wgs && n_smpls > 2000 then 0.001 else 0.0, # see https://odelaneau.github.io/shapeit5/docs/documentation/phase_rare/
          phase_extra_args = phase_extra_args,
          docker = docker_repository_with_sep + shapeit5_docker
      }
      if (length(batches) > 1 || defined(pedigree_file)) {
        call vcf_split {
          input:
            vcf_file = vcf_shapeit5.pgt_vcf_file,
            batches = batches,
            sample_id_file = sample_id_split_tsv.file,
            docker = docker_repository_with_sep + bcftools_docker
        }
      }
    }

    call vcf_concat {
      input:
        vcf_files = vcf_qc.xcl_vcf_file,
        filebase = sample_set_id + ".xcl",
        docker = docker_repository_with_sep + bcftools_docker
    }

    Array[Array[File]] batch_slices = if (length(batches) > 1 || defined(pedigree_file)) then transpose(select_all(vcf_split.vcf_files)) else [vcf_shapeit5.pgt_vcf_file]
    scatter (idx in range(n_batches)) {
      call vcf_ligate {
        input:
          vcf_files = batch_slices[idx],
          pedigree_file = if use_shapeit5_ligate then pedigree_file else None,
          filebase = pfxs[idx] + ".pgt",
          docker = docker_repository_with_sep + if use_shapeit5_ligate then shapeit5_docker else bcftools_docker
      }
      call vcf_import {
        input:
          pgt_vcf_file = vcf_ligate.vcf_file,
          pgt_vcf_idx = vcf_ligate.vcf_idx,
          unphased_vcf_file = if mode_is_vcf then data_paths_with_sep[idx] + batch_tbl["vcf"][idx] else select_first([unphased_vcf_files])[idx],
          unphased_vcf_idx = if mode_is_vcf then data_paths_with_sep[idx] + batch_tbl["vcf_index"][idx] else select_first([unphased_vcf_idxs])[idx],
          docker = docker_repository_with_sep + bcftools_docker
      }
    }
  }

  if (target == "calls" || target == "pngs") {
    scatter (idx in range(n_batches)) {
      call get_max_nrecords { input: vcf_idx = if mode == "pvcf" then data_paths_with_sep[idx] + batch_tbl["vcf_index"][idx] else select_first([vcf_import.vcf_idx])[idx], docker = docker_repository_with_sep + bcftools_docker }
      call vcf_mocha {
        input:
          n_smpls = if mode == "pvcf" then batch_tbl["n_smpls"][idx] else select_first([vcf_import.n_smpls])[idx],
          max_n_markers = get_max_nrecords.n,
          assembly = ref.name,
          rules_file = ref.rules_file,
          pvcf_file = if mode == "pvcf" then data_paths_with_sep[idx] + batch_tbl["vcf"][idx] else select_first([vcf_import.vcf_file])[idx],
          pvcf_idx = if mode == "pvcf" then data_paths_with_sep[idx] + batch_tbl["vcf_index"][idx] else select_first([vcf_import.vcf_idx])[idx],
          sample_tsv_file = select_first([sample_tsv.file, sample_tsv_file]),
          xcl_vcf_file = if mode == "pvcf" then data_paths_with_sep[idx] + batch_tbl["xcl_vcf"][idx] else vcf_concat.vcf_file,
          xcl_vcf_idx = if mode == "pvcf" then data_paths_with_sep[idx] + batch_tbl["xcl_vcf_index"][idx] else vcf_concat.vcf_idx,
          cnp_file = ref.cnp_file,
          mhc_reg = ref.mhc_reg,
          kir_reg = ref.kir_reg,
          mocha_extra_args = mocha_extra_args,
          ext_string = ext_string,
          docker = docker_repository_with_sep + bcftools_docker
      }
      if (target == "pngs") {
        call mocha_plot {
          input:
            vcf_file = vcf_mocha.mocha_vcf_file,
            vcf_idx = vcf_mocha.mocha_vcf_idx,
            stats_tsv = vcf_mocha.stats_tsv,
            calls_tsv = vcf_mocha.calls_tsv,
            cyto_file = ref.cyto_file,
            n_chrs = length(select_all(chrs)),
            call_rate_thr = if wgs then None else sample_call_rate_thr,
            baf_auto_thr = baf_auto_thr,
            mocha_plot_extra_args = mocha_plot_extra_args,
            wgs = wgs,
            docker = docker_repository_with_sep + r_mocha_docker
        }
      }
    }
    call tsv_concat as mocha_stats_tsv { input: tsv_files = vcf_mocha.stats_tsv, filebase = sample_set_id + ".stats", docker = basic_bash_docker }
    call tsv_concat as mocha_calls_tsv { input: tsv_files = vcf_mocha.calls_tsv, filebase = sample_set_id + ".calls", docker = basic_bash_docker }
    call mocha_summary {
      input:
        calls_tsv = mocha_calls_tsv.file,
        stats_tsv = mocha_stats_tsv.file,
        ucsc_beds = vcf_mocha.ucsc_bed,
        cyto_file = ref.cyto_file,
        n_chrs = length(select_all(chrs)),
        filebase = sample_set_id,
        call_rate_thr = if wgs then None else sample_call_rate_thr,
        baf_auto_thr = baf_auto_thr,
        docker = docker_repository_with_sep + r_mocha_docker
    }
  }

  # generate a table summarizing the main output files and serialize the table to disk
  # vcf_files and vcf_idxs are defined in the output section
  scatter (idx in range(n_batches)) {
    String basename_vcf_files = basename(vcf_files[idx])
    String basename_vcf_idxs = basename(vcf_idxs[idx])
    String basename_xcl_vcf_files = basename(if mode == "pvcf" then batch_tbl["xcl_vcf"][idx] else select_first([vcf_concat.vcf_file, ""]))
    String basename_xcl_vcf_idxs = basename(if mode == "pvcf" then batch_tbl["xcl_vcf_index"][idx] else select_first([vcf_concat.vcf_idx, ""]))
    String basename_pgt_vcf_files = basename(if mode != "pvcf" && target != "vcf" then select_first([vcf_ligate.vcf_file])[idx] else "")
    String basename_pgt_vcf_idxs = basename(if mode != "pvcf" && target != "vcf" then select_first([vcf_ligate.vcf_idx])[idx] else "")
  }
  Array[Pair[String,Array[String]]] output_pairs = flatten([[("batch_id", batch_tbl["batch_id"]),
    ("vcf", basename_vcf_files),
    ("vcf_index", basename_vcf_idxs)],
    if mode == "pvcf" || target == "vcf" then [] else
    [("n_smpls", select_first([vcf_import.n_smpls, unphased_vcf_n_smpls])),
    ("xcl_vcf", basename_xcl_vcf_files),
    ("xcl_vcf_index", basename_xcl_vcf_idxs),
    ("pgt_vcf", basename_pgt_vcf_files),
    ("pgt_vcf_index", basename_pgt_vcf_idxs)]])
  Map[String, Array[String]] output_map = as_map(output_pairs)
  # cannot use output_keys = keys(output_map) because of unresolved Cromwell bug
  # https://github.com/broadinstitute/cromwell/issues/5559
  scatter (p in output_pairs) {
    String output_keys = p.left
    Array[String] output_tsv_cols = output_map[p.left]
  }
  # this is run as a separate task rather than using write_tsv() as Cromwell can break the WDL specification
  # https://support.terra.bio/hc/en-us/community/posts/360071465631-write-lines-write-map-write-tsv-write-json-fail-when-run-in-a-workflow-rather-than-in-a-task
  call write_tsv {
    input:
      tsv = flatten([[output_keys], transpose(output_tsv_cols)]),
      filebase = sample_set_id + ".mocha",
      docker = basic_bash_docker
  }

  output {
    File? green_idat_tsv_file = green_idat_tsv.file
    File? red_idat_tsv_file = red_idat_tsv.file
    File? gtc_tsv_file = gtc_tsv.file
    File? cel_tsv_file = cel_tsv.file
    File? affy_tsv_file = affy_tsv.file
    File? stats_file = if defined(mocha_stats_tsv.file) || defined(sample_tsv.file) then select_first([mocha_stats_tsv.file, sample_tsv.file]) else None
    File? calls_file = mocha_calls_tsv.file
    File? ucsc_bed = mocha_summary.ucsc_bed
    File? summary_pdf = mocha_summary.summary_pdf
    File? pileup_pdf = mocha_summary.pileup_pdf
    Array[File]? png_files = if target == "pngs" then flatten(select_all(select_first([mocha_plot.png_files]))) else None
    Array[File]? gtc_files = if mode == "idat" && gtc_output then flatten(select_first([idat2gtc.gtc_files])) else None
    Array[File]? chp_files = if mode == "cel" && chp_output then flatten(select_first([cel2chp.chp_files])) else None
    Array[File]? snp_files = if mode == "cel" && chp_output then select_first([cel2chp.snp_file]) else None
    Array[File] vcf_files = select_first([vcf_mocha.mocha_vcf_file, vcf_import.vcf_file, unphased_vcf_files])
    Array[File] vcf_idxs = select_first([vcf_mocha.mocha_vcf_idx, vcf_import.vcf_idx, unphased_vcf_idxs])
    File? xcl_vcf_file = vcf_concat.vcf_file
    File? xcl_vcf_idx = vcf_concat.vcf_idx
    Array[File]? pgt_vcf_files = vcf_ligate.vcf_file
    Array[File]? pgt_vcf_idxs = vcf_ligate.vcf_idx
    File mocha_tsv_file = write_tsv.file
  }

  meta {
    author: "Giulio Genovese (with help from Chris Whelan)"
    email: "giulio.genovese@gmail.com"
    description: "See the [MoChA](https://github.com/freeseek/mocha) website for more information"
  }
}

task tsv_sorted {
  input {
    File tsv_file
    String column
    Boolean check_dups = false

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  String filebase = basename(tsv_file, ".tsv")

  command <<<
    set -euo pipefail
    mv "~{tsv_file}" .
    col=$(head -n1 "~{basename(tsv_file)}" | sed 's/\r$//' | tr '\t' '\n' | awk -F"\t" '$0=="~{column}" {print NR}')
    if [ "$col" == "" ]; then
      echo "Column \"~{column}\" does not exist" 1>&2
      exit 1
    fi
    cat "~{basename(tsv_file)}" | sed 's/\r$//' | (read -r; printf "%s\n" "$REPLY"; sort -k $col,$col -s -t $'\t' -T .) > "~{filebase}.sorted.tsv"
    ~{if check_dups then
      "dups=$(tail -n+2 \"" + filebase + ".sorted.tsv\" | sed 's/\\r$//' | cut -f$col | uniq --repeated)\n" +
      "if [[ ! -z \"$dups\" ]]; then\n  echo -e \"Duplicate " + column + "(s) found:\\n$dups\" 1>&2\n  exit 1\nfi"
    else ""}
    rm "~{basename(tsv_file)}"
  >>>

  output {
    File file = filebase + ".sorted.tsv"
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

task tsv_column {
  input {
    File tsv_file
    String column
    Boolean fail = true
    Boolean check_dups = false

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mv "~{tsv_file}" .
    col=$(head -n1 "~{basename(tsv_file)}" | sed 's/\r$//' | tr '\t' '\n' | awk -F"\t" '$0=="~{column}" {print NR}')
    if [ "$col" == "" ]; then
      ~{if fail then
        "echo \"Column \\\"" + column + "\\\" does not exist\" 1>&2\n" +
        "exit 1"
      else
        "touch \"" + column + ".lines\"\n" +
        "echo \"true\""}
    else
      ~{if check_dups then
        "dups=$(tail -n+2 \"" + basename(tsv_file) + "\" | sed 's/\\r$//' | cut -f$col | sort | uniq --repeated)\n" +
        "if [[ ! -z \"$dups\" ]]; then\n  echo -e \"Duplicate " + column + "(s) found:\\n$dups\" 1>&2\n  exit 1\nfi"
      else ""}
      ~{if column != "call_rate" then
        "tail -n+2 \"" + basename(tsv_file) + "\" | cut -f$col > \"" + column + ".lines\""
      else
        "max_call_rate=$(tail -n+2 \"" + basename(tsv_file) + "\" | sed 's/\\r$//' | cut -f$col | sort -g -T . | tail -n1)\n" +
        "tail -n+2 \"" + basename(tsv_file) + "\" | sed 's/\\r$//' | cut -f$col | if [[ $max_call_rate > 1.0 ]]; then awk '{print $0/100}'; else cat; fi > \"" + column + ".lines\"\n"}
      echo "false"
    fi
    rm "~{basename(tsv_file)}"
  >>>

  output {
    Boolean failed = read_boolean(stdout())
    File file = column + ".lines"
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

task tsv_concat {
  input {
    Array[File]+ tsv_files
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
    tsv_files=~{write_lines(tsv_files)}
    ~{if length(tsv_files) > 1 then
      "cat $tsv_files | tr '\\n' '\\0' | xargs -0 mv -t .\n" +
      "sed -i 's/^.*\\///' $tsv_files\n" +
      "(head -n1 \"" + basename(tsv_files[0]) + "\";\n" +
      "cat $tsv_files | tr '\\n' '\\0' | xargs -0 tail -qn+2) | sed 's/\\r$//'> \"" + filebase + ".tsv\"\n" +
      "cat $tsv_files | tr '\\n' '\\0' | xargs -0 rm"
    else "mv \"" + tsv_files[0] + "\" \"" + filebase + ".tsv\""}
  >>>

  output {
    File file = filebase + ".tsv"
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

task lst_flatten {
  input {
    Array[File]+ lst_files
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
    lst_files=~{write_lines(lst_files)}
    cat $lst_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $lst_files
    cat $lst_files | tr '\n' '\0' | xargs -0 awk 1 > "~{filebase}.lines"
    cat $lst_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File file = filebase + ".lines"
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

task lst_concat {
  input {
    Array[File]+ lst_files
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
    lst_files=~{write_lines(lst_files)}
    cat $lst_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $lst_files
    cat $lst_files | tr '\n' '\0' | xargs -0 -n 1 awk '{printf $0"\t"} END {printf "\n"}' | sed 's/\t$//' > "~{filebase}.tsv"
    cat $lst_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File file = filebase + ".tsv"
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

task lst_paste {
  input {
    Array[File]+ lst_files
    Array[String]+ headers
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
    headers=~{write_lines(headers)}
    lst_files=~{write_lines(lst_files)}
    cat $lst_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $lst_files
    (cat $headers | tr '\n' '\t' | sed 's/\t$/\n/'; cat $lst_files | tr '\n' '\0' | xargs -0 paste) > "~{filebase}.tsv"
    cat $lst_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File file = filebase + ".tsv"
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

# this task generates sub batches from batches and then returns them in the same order as batches
# batch_id_file should contain batch names without the delim string, or else behavior is undefined
task batch_scatter {
  input {
    File batch_id_file
    Int sub_batch_size
    String delim

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mv "~{batch_id_file}" .
    awk -F"\t" 'NR==FNR {x[$0]++} NR>FNR {n=x[$0]/int((x[$0]-1)/~{sub_batch_size}+1);
      if (x[$0]>n) print $0"~{delim}"int(y[$0]/n); else print $0; y[$0]++}' \
      "~{basename(batch_id_file)}" "~{basename(batch_id_file)}" > sub_batch_ids.lines
    uniq sub_batch_ids.lines
    uniq sub_batch_ids.lines | cut -d"~{delim}" -f1 | awk '!x[$0]++ {idx++} {print idx-1}' 1>&2
    rm "~{basename(batch_id_file)}"
  >>>

  output {
    File sub_batch_id = "sub_batch_ids.lines"
    Array[String] sub_batches = read_lines(stdout())
    Array[Int] idxs = read_lines(stderr())
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

task get_reheader_maps {
  input {
    Array[String]+ batches
    File barcodes_file
    File sample_id_file
    File batch_id_file

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    batches=~{write_lines(batches)}
    mv "~{barcodes_file}" .
    mv "~{sample_id_file}" .
    mv "~{batch_id_file}" .
    mkdir maps
    sed 's/^.*\///;s/.gz$//;s/_Grn.idat$//;s/.gtc$//;s/.CEL$//;s/.AxiomGT1.chp$//;s/.birdseed-v2.chp$//' "~{basename(barcodes_file)}" | \
      paste -d$'\t' - "~{basename(sample_id_file)}" | \
      paste -d$'\t' - "~{basename(batch_id_file)}" | \
      awk -F"\t" -v OFS="\t" '{print $1,$2>"maps/"$3".map"}'
    sed 's/^/maps\//;s/$/.map/' $batches
    rm "~{basename(barcodes_file)}"
    rm "~{basename(sample_id_file)}"
    rm "~{basename(batch_id_file)}"
  >>>

  output {
    Directory maps = "maps"
    # cannot use suffix(batches, ".map") because of unresolved Cromwell bug
    # https://github.com/broadinstitute/cromwell/issues/5549
    Array[File] map_files = read_lines(stdout())
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

task write_tsv {
  input {
    Array[Array[String]] tsv
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
    mv ~{write_tsv(tsv)} "~{filebase}.tsv"
  >>>

  output {
    File file = filebase + ".tsv"
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

task csv2bam {
  input {
    String plugin = "gtc2vcf"
    File csv_file
    File ref_fasta
    Array[File]+ ref_fasta_idxs

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float csv_size = (if basename(csv_file) != basename(csv_file, ".gz") then 4.0 else 1.0) * size(csv_file, "GiB")
  Float ref_size = size(ref_fasta, "GiB")
  Float index_size = size(ref_fasta_idxs, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * csv_size + ref_size + index_size)])
  # if gtc2vcf was memory efficient this requirement could be relaxed
  Float memory = select_first([memory_override, ceil(7.25 + 2.5 * csv_size)])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])

  command <<<
    set -euo pipefail
    echo "~{sep("\n", flatten([[csv_file, ref_fasta], ref_fasta_idxs]))}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    bcftools +~{plugin} \
      --csv "~{basename(csv_file)}" \
      --fasta-flank | \
      bwa mem~{if cpu > 1 then " -t " + cpu else ""} -M "~{basename(ref_fasta)}" - | \
      samtools view -bS -o "~{basename(basename(csv_file, ".gz"), ".csv")}.bam"
    echo "~{sep("\n", flatten([[csv_file, ref_fasta], ref_fasta_idxs]))}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File bam_file = basename(basename(csv_file, ".gz"), ".csv") + ".bam"
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

# https://support.terra.bio/hc/en-us/community/posts/360071476431-Terra-fails-to-delocalize-files-listed-through-read-lines-
task idat2gtc {
  input {
    File bpm_file
    File egt_file
    Array[File]+ green_idat_files
    Array[File]+ red_idat_files
    Int preset = 4
    String filebase

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0

    Float mult = 2.0 # to estimate the amount of memory required given the size of the manifest files
  }

  Float bpm_size = (if basename(bpm_file) != basename(bpm_file, ".gz") then 4.0 else 1.0) * size(bpm_file, "GiB")
  Float egt_size = (if basename(egt_file) != basename(egt_file, ".gz") then 2.0 else 1.0) * size(egt_file, "GiB")
  Float green_idat_size = length(green_idat_files) * size(green_idat_files[0], "GiB")
  Float red_idat_size = length(red_idat_files) * size(red_idat_files[0], "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + bpm_size + egt_size + 4.0 * (green_idat_size + red_idat_size))])
  Float memory = select_first([memory_override, 3.5 + mult * (bpm_size + egt_size)])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])

  command <<<
    set -euo pipefail
    green_idat_files=~{write_lines(green_idat_files)}
    red_idat_files=~{write_lines(red_idat_files)}
    mv "~{bpm_file}" .
    mv "~{egt_file}" .
    cat $green_idat_files $red_idat_files | tr '\n' '\0' | xargs -0 mv -t .
    ~{if basename(bpm_file) != basename(bpm_file, ".gz") then "gunzip --force \"" + basename(bpm_file) + "\"" else ""}
    ~{if basename(egt_file) != basename(egt_file, ".gz") then "gunzip --force \"" + basename(egt_file) + "\"" else ""}
    (grep -h "\.gz" $green_idat_files $red_idat_files || if [[ $? -eq 1 ]]; then true; else exit $?; fi) | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 gunzip --force
    sed -i 's/^.*\///;s/.gz$//' $green_idat_files
    sed -i 's/^.*\///;s/.gz$//' $red_idat_files
    bcftools +gtc2vcf --idat --gtcs $green_idat_files --output "~{filebase}.green_idat.tsv"
    bcftools +gtc2vcf --idat --gtcs $red_idat_files --output "~{filebase}.red_idat.tsv"
    bcftools +idat2gtc \
      --bpm "~{basename(bpm_file, ".gz")}" \
      --egt "~{basename(egt_file, ".gz")}" \
      --grn-idats $green_idat_files \
      --red-idats $red_idat_files \
      --output gtcs \
      --preset ~{preset} \
      --autocall-date ""
    sed 's/^/gtcs\//;s/_Grn\.idat$/.gtc/' $green_idat_files
    rm "~{basename(bpm_file, ".gz")}"
    rm "~{basename(egt_file, ".gz")}"
    cat $green_idat_files $red_idat_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File green_idat_tsv = filebase + ".green_idat.tsv"
    File red_idat_tsv = filebase + ".red_idat.tsv"
    Directory gtcs = "gtcs"
    Array[File] gtc_files = read_lines(stdout())
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

# https://support.terra.bio/hc/en-us/community/posts/360071476431-Terra-fails-to-delocalize-files-listed-through-read-lines-
task cel2affy {
  input {
    File? xml_file
    File? zip_file
    Array[File]+ cel_files
    File? cdf
    File? chrX_snps
    File? special_snps
    File? chrX_probes
    File? chrY_probes
    File? probeset_file
    Array[String]? chip_type
    Boolean table_output = false
    String? qmethod_spec
    File? read_models_brlmmp
    String? set_analysis_name
    File? target_sketch
    String? set_gender_method
    String? em_gender
    Float? female_thresh
    Float? male_thresh
    String filebase

    String docker
    Int? cpu_override # unfortunately this task is not parallelizable
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 0
    Int maxRetries = 0
  }

  Float xml_size = size(xml_file, "GiB")
  Float zip_size = size(zip_file, "GiB")
  Float cel_size = length(cel_files) * size(cel_files[0], "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + xml_size + zip_size + (if table_output then 8.0 else 6.0) * cel_size)])
  Float memory = select_first([memory_override, 2 * 7.25])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])

  String? zip_dir = if defined(zip_file) then basename(select_first([zip_file]), ".zip") else None
  String analysis_name = if defined(set_analysis_name) then set_analysis_name else "AxiomGT1"

  command <<<
    set -euo pipefail
    cel_files=~{write_lines(cel_files)}
    echo "~{sep("\n", select_all([xml_file, zip_file, cdf, chrX_snps, special_snps, chrX_probes, chrY_probes, read_models_brlmmp, target_sketch]))}" | \
      cat - $cel_files | tr '\n' '\0' | xargs -0 mv -t .
    (grep "\.gz" $cel_files  || if [[ $? -eq 1 ]]; then true; else exit $?; fi) | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 gunzip --force
    sed -i 's/^.*\///;s/.gz$//' $cel_files
    ~{if defined(zip_file) then "unzip -jd \"" + zip_dir + "\" \"" + basename(select_first([zip_file])) + "\" 1>&2" else ""}
    bcftools +affy2vcf --cel --chps $cel_files --output "~{filebase}.cel.tsv"
    echo "cel_files" | cat - $cel_files > cel_files.lines
    apt-probeset-genotype \
      ~{if defined(zip_dir) then "--analysis-files-path \"" + zip_dir + "\"" else ""} \
      ~{if defined(xml_file) then "--xml-file \"" + basename(select_first([xml_file])) + "\"" else ""} \
      --cel-files cel_files.lines \
      ~{if defined(cdf) then "--cdf-file \"" + basename(select_first([cdf])) + "\"" else ""} \
      ~{if defined(chrX_snps) then "--chrX-snps \"" + basename(select_first([chrX_snps])) + "\"" else ""} \
      ~{if defined(special_snps) then "--special-snps \"" + basename(select_first([special_snps])) + "\"" else ""} \
      ~{if defined(chrX_probes) then "--chrX-probes \"" + basename(select_first([chrX_probes])) + "\"" else ""} \
      ~{if defined(chrY_probes) then "--chrY-probes \"" + basename(select_first([chrY_probes])) + "\"" else ""} \
      ~{if defined(probeset_file) then "--probeset-ids \"" + select_first([probeset_file]) + "\"" else ""} \
      ~{if defined(chip_type) then "--chip-type " + sep(" --chip-type ", select_first([chip_type])) else ""} \
      --table-output ~{table_output} \
      --cc-chp-output \
      ~{if table_output then "--summaries" else ""} \
      ~{if defined(qmethod_spec) then "--qmethod-spec " + qmethod_spec else ""} \
      ~{if defined(read_models_brlmmp) then "--read-models-brlmmp \"" + basename(select_first([read_models_brlmmp])) + "\"" else ""} \
      --write-models \
      ~{if defined(target_sketch) then "--target-sketch \"" + basename(select_first([target_sketch])) + "\"" else ""} \
      ~{if defined(set_analysis_name) then "--set-analysis-name " + set_analysis_name else ""} \
      ~{if defined(set_gender_method) then "--set-gender-method " + set_gender_method else ""} \
      ~{if defined(em_gender) then "--em-gender " + em_gender else ""} \
      ~{if defined(female_thresh) then "--female-thresh " + female_thresh else ""} \
      ~{if defined(male_thresh) then "--male-thresh " + male_thresh else ""}
    for sfx in snp-posteriors report~{if table_output then " calls confidences summary normalized-summary" else ""}; do
      mv "~{analysis_name}.$sfx.txt" "~{filebase}.$sfx.txt"
    done
    sed 's/^/cc-chp\//;s/\.CEL$/.~{analysis_name}.chp/' $cel_files
    echo "~{sep("\n", select_all([xml_file, zip_file, cdf, chrX_snps, special_snps, chrX_probes, chrY_probes, read_models_brlmmp, target_sketch]))}" | \
      sed 's/^.*\///' | cat - $cel_files | tr '\n' '\0' | xargs -0 rm
    ~{if defined(zip_file) then "rm -r \"" + zip_dir + "\"" else ""}
    rm cel_files.lines
  >>>

  output {
    File cel_tsv = filebase + ".cel.tsv"
    Directory chps = "cc-chp"
    Array[File] chp_files = read_lines(stdout())
    File log_file = "apt-probeset-genotype.log"
    File snp_file = filebase + ".snp-posteriors.txt"
    File report_file = filebase + ".report.txt"
    File? calls_file = if table_output then filebase + ".calls.txt" else None
    File? confidences_file = if table_output then filebase + ".confidences.txt" else None
    File? summary_file = if table_output then filebase + ".summary.txt" else None
    File? normalized_summary_file = if table_output then filebase + ".normalized-summary.txt" else None
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

task gtc2vcf {
  input {
    String? tags
    File bpm_file
    File csv_file
    File egt_file
    File ref_fasta
    File ref_fasta_fai
    Int? gc_window_size
    Array[File]+ gtc_files
    Int capacity = 32768
    String? gtc2vcf_extra_args
    File? sam_file
    File? reheader_file
    String filebase

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float bpm_size = (if basename(bpm_file) != basename(bpm_file, ".gz") then 4.0 else 1.0) * size(bpm_file, "GiB")
  Float csv_size = (if basename(csv_file) != basename(csv_file, ".gz") then 4.0 else 1.0) * size(csv_file, "GiB")
  Float egt_size = (if basename(egt_file) != basename(egt_file, ".gz") then 2.0 else 1.0) * size(egt_file, "GiB")
  Float ref_size = size(ref_fasta, "GiB")
  Float gtc_size = length(gtc_files) * size(gtc_files[0], "GiB")
  Float sam_size = size(sam_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + bpm_size + csv_size + egt_size + ref_size + 8.0 * gtc_size + sam_size)])
  # due to heavy random access to the reference genome, it is important here that enough memory to cache the reference is provided
  Float memory = select_first([memory_override, 3.5 + 2.0 * csv_size + 2.0 * egt_size + ref_size + 19.0 * length(gtc_files) * capacity / 1024 / 1024 / 1024])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])

  command <<<
    set -euo pipefail
    gtc_files=~{write_lines(gtc_files)}
    echo "~{sep("\n", select_all([bpm_file, csv_file, egt_file, ref_fasta, ref_fasta_fai, sam_file, reheader_file]))}" | \
      cat - $gtc_files | tr '\n' '\0' | xargs -0 mv -t .
    ~{if basename(bpm_file) != basename(bpm_file, ".gz") then "gunzip --force \"" + basename(bpm_file) + "\"" else ""}
    ~{if basename(egt_file) != basename(egt_file, ".gz") then "gunzip --force \"" + basename(egt_file) + "\"" else ""}
    (grep "\.gz" $gtc_files || if [[ $? -eq 1 ]]; then true; else exit $?; fi) | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 gunzip --force
    sed -i 's/^.*\///;s/.gz$//' $gtc_files
    ~{if defined(reheader_file) then "sed -i 's/ /\\\\ /g' \"" + basename(select_first([reheader_file])) + "\"" else ""}
    bcftools +gtc2vcf \
      --no-version \
      --output-type u \
      ~{if defined(tags) then "--tags " + tags else ""} \
      --bpm "~{basename(bpm_file, ".gz")}" \
      --csv "~{basename(csv_file)}" \
      --egt "~{basename(egt_file, ".gz")}" \
      --fasta-ref "~{basename(ref_fasta)}" \
      ~{if defined(gc_window_size) then "--gc-window-size " + select_first([gc_window_size]) else ""} \
      --gtcs $gtc_files \
      ~{if capacity != 32768 then "--capacity " + capacity else ""} \
      ~{if defined(gtc2vcf_extra_args) then select_first([gtc2vcf_extra_args]) else ""} \
      --extra "~{filebase}.gtc.tsv" \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      ~{if defined(sam_file) then "--sam-flank \"" + basename(select_first([sam_file])) + "\"" else ""} | \
      ~{if defined(reheader_file) then "bcftools reheader --samples \"" + basename(select_first([reheader_file])) + "\" |" else ""} \
      bcftools sort --output-type u --temp-dir ./bcftools. | \
    bcftools norm \
      --no-version \
      --output "~{filebase}.bcf" \
      --output-type b \
      --check-ref x \
      --fasta-ref "~{basename(ref_fasta)}" \
      ~{if cpu > 1 then " --threads " + (cpu - 1) else ""} \
      --write-index
    bcftools query --list-samples "~{filebase}.bcf" | tee "~{filebase}.sample_id.lines" | wc -l
    rm "~{basename(bpm_file, ".gz")}"
    rm "~{basename(egt_file, ".gz")}"
    echo "~{sep("\n", select_all([csv_file, ref_fasta, ref_fasta_fai, sam_file, reheader_file]))}" | \
      sed 's/^.*\///' | cat - $gtc_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File gtc_tsv = filebase + ".gtc.tsv"
    File vcf_file = filebase + ".bcf"
    File vcf_idx = filebase + ".bcf.csi"
    File sample_id_lines = filebase + ".sample_id.lines"
    Int n_smpls = read_int(stdout())
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

task chp2vcf {
  input {
    String? tags
    File csv_file
    File ref_fasta
    File ref_fasta_fai
    Int? gc_window_size
    File? probeset_file
    File snp_file
    Array[File]+ chp_files
    File? sam_file
    File? reheader_file
    String filebase

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float csv_size = (if basename(csv_file) != basename(csv_file, ".gz") then 4.0 else 1.0) * size(csv_file, "GiB")
  Float ref_size = size(ref_fasta, "GiB")
  Float snp_size = (if basename(snp_file) != basename(snp_file, ".gz") then 2.0 else 1.0) * size(snp_file, "GiB")
  Float chp_size = length(chp_files) * size(chp_files[0], "GiB")

  Float sam_size = size(sam_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + csv_size + ref_size + 8.0 * chp_size + snp_size + sam_size)])
  # due to heavy random access to the reference genome, it is important here that enough memory to cache the reference is provided
  Float memory = select_first([memory_override, 3.5 + 2.0 * (csv_size + snp_size) + ref_size + length(chp_files) * 32768 / 1024 / 1024 / 1024])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])

  command <<<
    set -euo pipefail
    chp_files=~{write_lines(chp_files)}
    echo "~{sep("\n", select_all([csv_file, ref_fasta, ref_fasta_fai, snp_file, sam_file, reheader_file]))}" | \
      cat - $chp_files | tr '\n' '\0' | xargs -0 mv -t .
    (grep "\.gz" $chp_files || if [[ $? -eq 1 ]]; then true; else exit $?; fi) | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 gunzip --force
    sed -i 's/^.*\///;s/.gz$//' $chp_files
    ~{if defined(reheader_file) then "sed -i 's/ /\\\\ /g' \"" + basename(select_first([reheader_file])) + "\"" else ""}
    bcftools +affy2vcf \
      --no-version \
      --output-type u \
      ~{if defined(tags) then "--tags " + tags else ""} \
      --csv "~{basename(csv_file)}" \
      --fasta-ref "~{basename(ref_fasta)}" \
      ~{if defined(gc_window_size) then "--gc-window-size " + select_first([gc_window_size]) else ""} \
      ~{if defined(probeset_file) then "--probeset-ids \"" + select_first([probeset_file]) + "\"" else ""} \
      --snp "~{basename(snp_file)}" \
      --chps $chp_files \
      --extra "~{filebase}.affy.tsv" \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      ~{if defined(sam_file) then "--sam-flank \"" + basename(select_first([sam_file])) + "\"" else ""} | \
      ~{if defined(reheader_file) then "bcftools reheader --samples \"" + basename(select_first([reheader_file])) + "\" |" else ""} \
    bcftools sort --output-type u --temp-dir ./bcftools. | \
    bcftools norm \
      --no-version \
      --output "~{filebase}.bcf" \
      --output-type b \
      --check-ref x \
      --fasta-ref "~{basename(ref_fasta)}" \
      ~{if cpu > 1 then " --threads " + (cpu - 1) else ""} \
      --write-index
    bcftools query --list-samples "~{filebase}.bcf" | tee "~{filebase}.sample_id.lines" | wc -l
    echo "~{sep("\n", select_all([csv_file, ref_fasta, ref_fasta_fai, snp_file, sam_file, reheader_file]))}" | \
      sed 's/^.*\///' | cat - $chp_files | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File affy_tsv = filebase + ".affy.tsv"
    File vcf_file = filebase + ".bcf"
    File vcf_idx = filebase + ".bcf.csi"
    File sample_id_lines = filebase + ".sample_id.lines"
    Int n_smpls = read_int(stdout())
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

task txt2vcf {
  input {
    String? tags
    File csv_file
    File ref_fasta
    File ref_fasta_fai
    Int? gc_window_size
    File? probeset_file
    File calls_file
    File? confidences_file
    File summary_file
    File report_file
    File snp_file
    File? sam_file
    File? reheader_file
    String filebase

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float csv_size = (if basename(csv_file) != basename(csv_file, ".gz") then 4.0 else 1.0) * size(csv_file, "GiB")
  Float ref_size = size(ref_fasta, "GiB")
  Float calls_size = (if basename(calls_file) != basename(calls_file, ".gz") then 2.0 else 1.0) * size(calls_file, "GiB")
  Float confidences_size = (if defined(confidences_file) && basename(select_first([confidences_file])) != basename(select_first([confidences_file]), ".gz") then 2.0 else 1.0) * size(confidences_file, "GiB")
  Float summary_size = (if basename(summary_file) != basename(summary_file, ".gz") then 2.0 else 1.0) * size(summary_file, "GiB")
  Float report_size = (if basename(report_file) != basename(report_file, ".gz") then 2.0 else 1.0) * size(report_file, "GiB")
  Float snp_size = (if basename(snp_file) != basename(snp_file, ".gz") then 2.0 else 1.0) * size(snp_file, "GiB")

  Float sam_size = size(sam_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + csv_size + ref_size + 8.0 * (calls_size + confidences_size + summary_size) + 2.0 * report_size + snp_size + sam_size)])
  # due to heavy random access to the reference genome, it is important here that enough memory to cache the reference is provided
  Float memory = select_first([memory_override, 3.5 + 2.0 * csv_size + 2.0 * snp_size + ref_size])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])

  command <<<
    set -euo pipefail
        echo "~{sep("\n", select_all([csv_file, ref_fasta, ref_fasta_fai, calls_file, confidences_file, summary_file, report_file, snp_file, sam_file, reheader_file]))}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    ~{if basename(report_file) != basename(report_file, ".gz") then "gunzip --force \"" + basename(report_file) + "\"" else ""}
    (~{if basename(calls_file) != basename(calls_file, ".gz") then "z" else ""}grep -v ^# "~{basename(calls_file)}" || if [[ $? -eq 141 ]]; then true; else exit $?; fi) | \
      head -n1 | tr '\t' '\n' | tail -n+2 | \
      awk -F"\t" 'NR==FNR && $0!~"^#" {if ($1=="cel_files") print; else x[$1]=$0} NR>FNR {print x[$1]}' \
      "~{basename(report_file, ".gz")}" - > ~{filebase}.affy.tsv
    rm "~{basename(report_file, ".gz")}"
    ~{if defined(reheader_file) then "sed -i 's/ /\\\\ /g' \"" + basename(select_first([reheader_file])) + "\"" else ""}
    bcftools +affy2vcf \
      --no-version \
      --output-type u \
      ~{if defined(tags) then "--tags " + tags else ""} \
      --csv "~{basename(csv_file)}" \
      --fasta-ref "~{basename(ref_fasta)}" \
      ~{if defined(gc_window_size) then "--gc-window-size " + select_first([gc_window_size]) else ""} \
      ~{if defined(probeset_file) then "--probeset-ids \"" + select_first([probeset_file]) + "\"" else ""} \
      --calls "~{basename(calls_file)}" \
      ~{if defined(confidences_file) then "--confidences \"" + basename(select_first([confidences_file])) + "\"" else ""} \
      --summary "~{basename(summary_file)}" \
      --snp "~{basename(snp_file)}" \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      ~{if defined(sam_file) then "--sam-flank \"" + basename(select_first([sam_file])) + "\"" else ""} | \
    ~{if defined(reheader_file) then "bcftools reheader --samples \"" + basename(select_first([reheader_file])) + "\" |" else ""} \
    bcftools sort --output-type u --temp-dir ./bcftools. | \
    bcftools norm \
      --no-version \
      --output "~{filebase}.bcf" \
      --output-type b \
      --check-ref x \
      --fasta-ref "~{basename(ref_fasta)}" \
      ~{if cpu > 1 then " --threads " + (cpu - 1) else ""} \
      --write-index
    bcftools query --list-samples "~{filebase}.bcf" | tee "~{filebase}.sample_id.lines" | wc -l
    echo "~{sep("\n", select_all([csv_file, ref_fasta, ref_fasta_fai, calls_file, confidences_file, summary_file, snp_file, sam_file, reheader_file]))}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File affy_tsv = filebase + ".affy.tsv"
    File vcf_file = filebase + ".bcf"
    File vcf_idx = filebase + ".bcf.csi"
    File sample_id_lines = filebase + ".sample_id.lines"
    Int n_smpls = read_int(stdout())
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
        cm_mids = (win_size - 0.5 * ~{overlap_size_cm}) * np.arange(1, n_win - 1)
        pos_begs = np.concatenate(([0], 0 + np.interp(cm_begs, df_group['CM'], df_group['POS'], period = np.inf).astype(int)))
        pos_ends = np.concatenate((np.interp(cm_ends, df_group['CM'], df_group['POS'], period = np.inf).astype(int), [chr2len[fai_chr]]))
        pos_begs2 = np.concatenate(([0], 0 + np.interp(cm_begs + 0.5 * ~{overlap_size_cm}, df_group['CM'], df_group['POS'], period = np.inf).astype(int)))
        pos_ends2 = np.concatenate((np.interp(cm_ends - 0.5 * ~{overlap_size_cm}, df_group['CM'], df_group['POS'], period = np.inf).astype(int), [chr2len[fai_chr]]))
        df_out[fai_chr] = pd.DataFrame.from_dict({'CHR': fai_chr, 'BEG': pos_begs, 'END': pos_ends, 'BEG2': pos_begs2, 'END2': pos_ends2})
    df = pd.concat([df_out[fai_chr] for fai_chr in chr2len.keys()])
    df[['CHR', 'BEG', 'END', 'BEG2', 'END2']].to_csv('ref_scatter.tsv', sep='\t', header = False, index = False)
    CODE
    rm chr2len.tsv
    rm "~{basename(genetic_map_file)}"
  >>>

  output {
    File intervals_tsv = "ref_scatter.tsv"
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

# https://support.terra.bio/hc/en-us/community/posts/360071476431-Terra-fails-to-delocalize-files-listed-through-read-lines-
task vcf_scatter {
  input {
    File vcf_file
    File intervals_tsv # zero-based intervals
    Int clevel = 2

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
    mv "~{vcf_file}" .
    mv "~{intervals_tsv}" .
    awk -F"\t" '{print $1":"1+$2"-"$3"\t"NR-1}' "~{basename(intervals_tsv)}" > regions.lines
    bcftools query --list-samples "~{basename(vcf_file)}" | tee "~{filebase}.sample_id.lines" | wc -l > n_smpls.int
    bcftools annotate \
      --no-version \
      --output-type u \
      --remove ID,QUAL,INFO,^FMT/GT \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      "~{basename(vcf_file)}" | \
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
    rm "~{basename(vcf_file)}"
    rm "~{basename(intervals_tsv)}"
    rm regions.lines
  >>>

  output {
    File sample_id_lines = filebase + ".sample_id.lines"
    Int n_smpls = read_int("n_smpls.int")
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
    String filebase
    Boolean wgs = false

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
      "  --output \"" + filebase + ".bcf\" \\\n" +
      "  --output-type b \\\n" +
      (if wgs then "--missing-to-ref \\\n" else "") +
      "  --file-list $vcf_files \\\n" +
      "  --merge none \\\n" +
      "  --no-index \\\n" +
      (if cpu > 1 then "  --threads " + (cpu - 1) else "") + " \\\n" +
      "--write-index"
    else "mv \"" + vcf_files[0] + "\" \"" + filebase + ".bcf\"\n" +
      "bcftools index --force \"" + filebase + ".bcf\""}
    bcftools query --list-samples "~{filebase}.bcf" | wc -l
    ~{if length(vcf_files) > 1 then "cat $vcf_files | tr '\\n' '\\0' | xargs -0 rm" else ""}
  >>>

  output {
    File vcf_file = filebase + ".bcf"
    File vcf_idx = filebase + ".bcf.csi"
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

task vcf_qc {
  input {
    File vcf_file
    File? vcf_idx
    File? dup_file
    File sample_tsv_file
    Float? sample_call_rate_thr
    Float variant_call_rate_thr = 0.97
    Float dup_divergence_thr = 0.02
    Float genotype_exc_het_thr = 0.000001
    File? duplicate_samples_file
    File? extra_xcl_vcf_file

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + vcf_size)])
  Float memory = select_first([memory_override, 3.5])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])
  String filebase = basename(basename(vcf_file, ".bcf"), ".vcf.gz")

  command <<<
    set -euo pipefail
    echo "~{sep("\n", select_all([vcf_file, vcf_idx, dup_file, sample_tsv_file, duplicate_samples_file, extra_xcl_vcf_file]))}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    awk -F"\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
      NR>1 ~{if defined(sample_call_rate_thr) then "&& $(f[\"call_rate\"])<" + select_first([sample_call_rate_thr]) else ""} {print $(f["sample_id"])}' \
      "~{basename(sample_tsv_file)}" ~{if defined(duplicate_samples_file) then "| \\\n" +
      "  cat - \"" + basename(select_first([duplicate_samples_file])) + "\" | \\\n" +
      "  sort -T . | uniq " else ""}> samples_xcl.lines
    ~{if defined(vcf_idx) then "" else "bcftools index --force \"" + basename(vcf_file) + "\""}
    ~{if defined(dup_file) then
    "echo '##INFO=<ID=JK,Number=1,Type=Float,Description=\"Jukes Cantor\">' | \\\n" +
    "bcftools annotate --no-version --output-type u --annotations \"" + basename(select_first([dup_file])) + "\" --columns CHROM,FROM,TO,JK --header-lines /dev/stdin" + (if cpu > 1 then " --threads " + (cpu - 1) else "") + " \"" + basename(vcf_file) + "\" | \\"
    else ""}
    bcftools view --no-version --output-type u~{if cpu > 1 then " --threads " + (cpu - 1) else ""} --samples-file ^samples_xcl.lines --force-samples ~{if defined(dup_file) then "" else "\"" + basename(vcf_file) + "\""} | \
    bcftools +fill-tags --no-version --output-type u --targets ^Y,MT,chrY,chrM~{if cpu > 1 then " --threads " + (cpu - 1) else ""} -- --tags ExcHet,F_MISSING | \
    bcftools view --no-version --output-type u --drop-genotypes | \
    bcftools annotate --no-version --output "~{filebase}.xcl.bcf" --output-type b \
      --include 'FILTER!="." && FILTER!="PASS" ~{if defined(dup_file) then "|| INFO/JK<" + dup_divergence_thr else ""}|| INFO/ExcHet<~{genotype_exc_het_thr} || INFO/F_MISSING>1-~{variant_call_rate_thr}' \
      --remove ~{if defined(dup_file) then "^INFO/JK," else ""}^INFO/ExcHet,^INFO/F_MISSING~{if cpu > 1 then " --threads " + (cpu - 1) else ""} --write-index
    ~{if defined(extra_xcl_vcf_file) then "mv \"" + filebase + ".xcl.bcf\" \"" + filebase + ".tmp.bcf\"\n" +
      "mv \"" + filebase + ".xcl.bcf.csi\" \"" + filebase + ".tmp.bcf.csi\"\n" +
      "bcftools index --force \"" + basename(select_first([extra_xcl_vcf_file])) + "\"\n" +
      "bcftools merge --no-version --output \"" + filebase + ".xcl.bcf\" --output-type b --merge none \"" + filebase + ".tmp.bcf\" \"" + basename(select_first([extra_xcl_vcf_file])) + "\" --write-index\n" +
      "rm \"" + filebase + ".tmp.bcf\" \"" + filebase + ".tmp.bcf.csi\"" else ""}
    bcftools isec \
      --no-version \
      --output-type u \
      --complement \
      --exclude "N_ALT>1" \
      --write 1 \
      "~{basename(vcf_file)}" \
      "~{filebase}.xcl.bcf" | \
    bcftools view --no-version --output "~{filebase}.qc.bcf" --output-type b --min-ac 0 --exclude-uncalled --write-index
    bcftools index --nrecords "~{filebase}.qc.bcf"
    ~{if defined(vcf_idx) then "" else "rm \"" + basename(vcf_file) + ".csi\""}
    echo "~{sep("\n", select_all([vcf_file, vcf_idx, dup_file, sample_tsv_file, duplicate_samples_file, extra_xcl_vcf_file]))}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
    rm samples_xcl.lines
  >>>

  output {
    File xcl_vcf_file = filebase + ".xcl.bcf"
    File xcl_vcf_idx = filebase + ".xcl.bcf.csi"
    File qc_vcf_file = filebase + ".qc.bcf"
    File qc_vcf_idx = filebase + ".qc.bcf.csi"
    Int n_markers = read_int(stdout())
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

# hack https://github.com/samtools/bcftools/issues/1425 is employed in the end to fix the header
# the command requires BCFtools 1.14 due to bug https://github.com/samtools/bcftools/issues/1497
task vcf_shapeit5 {
  input {
    Int n_smpls
    Int n_markers
    File unphased_vcf_file
    File unphased_vcf_idx
    File genetic_map_file
    Int n_chrs
    File? pedigree_file
    File? sample_tsv_file
    Int? n_panel_smpls
    File? ref_vcf_file
    File? ref_vcf_idx
    File? ref_fasta_fai
    String input_region
    String scaffold_region
    String chr
    Boolean mhc = false # requires additional resources if phasing the MHC region
    Float maf = 0.0
    String? phase_extra_args

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? mult_override
    Float? memory_override
    Int? preemptible_override
    Int maxRetries = 0
  }

  Float vcf_size = size(unphased_vcf_file, "GiB")
  Float ref_size = size(ref_vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 3.0 * vcf_size + ref_size)])
  Float mult = select_first([mult_override, 12.0 * (if chr == "X" || chr == "chrX" then 1.5 else 1.0)])
  Float memory = select_first([memory_override, 3.5 + mult * n_markers * (n_smpls + select_first([n_panel_smpls, 0])) / 1024 / 1024 / 1024])
  Int cpu = select_first([cpu_override, 2 * ceil(memory / 13)]) * (if mhc then 2 else 1) # always require at least two CPUs (four for MHC windows)
  Int preemptible = select_first([preemptible_override, if mhc then 0 else 1]) # as the MHC phases slowly, do not run as preemptible

  String filebase = basename(basename(unphased_vcf_file, ".bcf"), ".vcf.gz")
  String dollar = "$"

  command <<<
    set -euo pipefail
    echo "~{sep("\n", select_all([unphased_vcf_file, unphased_vcf_idx, genetic_map_file, pedigree_file, sample_tsv_file, ref_vcf_file, ref_vcf_idx]))}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    ~{if defined(ref_fasta_fai) then "mv \"" + select_first([ref_fasta_fai]) + "\" ." else ""}
    chr=~{chr}; zcat "~{basename(genetic_map_file)}" | \
      sed 's/^~{n_chrs}/X/' | awk -v chr=${chr#chr} '$1==chr {print $2,$3,$4}' > genetic_map.txt
    ~{if defined(sample_tsv_file) then
    "awk -F\"\t\" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}\n" +
    "  NR>1 && $(f[\"computed_gender\"])==\"M\" {print $(f[\"sample_id\"])}' \\\n" +
    "  > \"" + filebase + ".haploids.lines\" \"" + basename(select_first([sample_tsv_file])) + "\""
    else ""}
    phase_common \
      ~{if cpu > 1 then "--thread " + cpu else ""} \
      --input "~{basename(unphased_vcf_file)}" \
      ~{if defined(ref_vcf_file) then "--reference \"" + basename(select_first([ref_vcf_file])) + "\"" else ""} \
      --map genetic_map.txt \
      ~{if defined(pedigree_file) then "--pedigree \"" + basename(select_first([pedigree_file])) + "\"" else ""} \
      ~{if defined(sample_tsv_file) then "--haploids \"" + filebase + ".haploids.lines\"" else ""} \
      --region ~{scaffold_region} \
      ~{if maf > 0 then "--filter-maf " + maf else ""} \
      --output "~{filebase}.~{if maf > 0 then "scaffold" else "pgt"}.bcf" \
      ~{if defined(phase_extra_args) then phase_extra_args else ""} \
      1>&2
    ~{if maf > 0 then
      "phase_rare \\\n" +
      (if cpu > 1 then "  --thread " + cpu + " \\\n" else "") +
      "  --input \"" + basename(unphased_vcf_file) + "\" \\\n" +
      "  --scaffold \"" + filebase + ".scaffold.bcf\" \\\n" +
      "  --map genetic_map.txt \\\n" +
      (if defined(pedigree_file) then "  --pedigree \"" + basename(select_first([pedigree_file])) + "\" \\\n" else "") +
      (if defined(sample_tsv_file) then "--haploids \"" + filebase + ".haploids.lines\" \\\n" else "") +
      "  --input-region " + input_region + " \\\n" +
      "  --scaffold-region " + scaffold_region + " \\\n" +
      "  --output \"" + filebase + ".pgt.bcf\" \\\n" +
      "1>&2\n" +
      "rm \"" + filebase + ".scaffold.bcf\"\n"
    else ""}rm genetic_map.txt~{if defined(sample_tsv_file) then " \"" + filebase + ".haploids.lines\"" else ""}
    ~{if defined(ref_fasta_fai) then
      "(echo -en \"##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t\"\n" +
      "bcftools query -l \"" + filebase + ".pgt.bcf\" | tr '\\n' '\\t' | sed 's/\\t$/\\n/') > tmp.vcf\n" +
      "bcftools reheader --fai \"" + basename(select_first([ref_fasta_fai])) + "\" --output fai.vcf --temp-prefix ./bcftools. tmp.vcf\n" +
      "mv \"" + filebase + ".pgt.bcf\" \"" + filebase + ".tmp.bcf\"\n" +
      "bcftools concat \\\n" +
      "  --no-version \\\n" +
      "  --output-type b \\\n" +
      "  --output \"" + filebase + ".pgt.bcf\" \\\n" +
      "  fai.vcf \"" + filebase + ".tmp.bcf\"\n" +
      "rm \"" + filebase + ".tmp.bcf\" fai.vcf tmp.vcf \"" + basename(select_first([ref_fasta_fai])) + "\""
      else ""}
    bcftools index --force "~{filebase}.pgt.bcf"
    echo "~{sep("\n", select_all([unphased_vcf_file, unphased_vcf_idx, genetic_map_file, pedigree_file, sample_tsv_file, ref_vcf_file, ref_vcf_idx]))}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File pgt_vcf_file = filebase + ".pgt.bcf"
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

# https://support.terra.bio/hc/en-us/community/posts/360071476431-Terra-fails-to-delocalize-files-listed-through-read-lines-
task vcf_split {
  input {
    File vcf_file
    Array[String]+ batches
    File sample_id_file
    Int clevel = 2

    String docker
    Int? disk_size_override
    Int? cpu_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])
  String filebase = basename(basename(vcf_file, ".bcf"), ".vcf.gz")

  command <<<
    set -euo pipefail
    filebases=~{write_lines(prefix(filebase + '.', batches))}
    mv "~{vcf_file}" .
    mv "~{sample_id_file}" .
    mkdir mendel
    sed 's/ /\\ /g;s/\t/,/g;s/$/\t-/' "~{basename(sample_id_file)}" | paste -d $'\t' - $filebases > samples_file.txt
    bcftools +split --keep-tags FMT/GT --output-type b~{clevel} --output vcfs --samples-file samples_file.txt "~{basename(vcf_file)}"
    cut -f1 $filebases | sed 's/^/vcfs\//;s/$/.bcf/'
    rm "~{basename(vcf_file)}"
    rm "~{basename(sample_id_file)}"
    rm samples_file.txt
  >>>

  output {
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

task vcf_concat {
  input {
    Array[File]+ vcf_files
    String filebase

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float vcf_size = size(vcf_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])

  command <<<
    set -euo pipefail
    vcf_files=~{write_lines(vcf_files)}
    cat $vcf_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $vcf_files
    cat $vcf_files | tr '\n' '\0' | xargs -0 -n 1 bcftools index --force
    bcftools concat \
      --no-version \
      --output "~{filebase}.bcf" \
      --output-type b \
      --allow-overlaps \
      --file-list $vcf_files \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      --write-index
    cat $vcf_files | tr '\n' '\0' | xargs -0 rm
    cat $vcf_files | sed 's/$/.csi/' | tr '\n' '\0' | xargs -0 rm
    bcftools query --list-samples "~{filebase}.bcf" | wc -l
  >>>

  output {
    File vcf_file = filebase + ".bcf"
    File vcf_idx = filebase + ".bcf.csi"
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


task vcf_ligate {
  input {
    Array[File]+ vcf_files
    File? pedigree_file
    String filebase
    Boolean use_shapeit5_ligate = false

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0

    Float mult = 12.0
  }

  Float vcf_size = size(vcf_files, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5 + (if use_shapeit5_ligate then 2.0 * mult else mult) * vcf_size / length(vcf_files)])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])

  command <<<
    set -euo pipefail
    ~{if defined(pedigree_file) then "mv \"" + select_first([pedigree_file]) + "\" ." else ""}
    vcf_files=~{write_lines(vcf_files)}
    cat $vcf_files | tr '\n' '\0' | xargs -0 mv -t .
    sed -i 's/^.*\///' $vcf_files
    cat $vcf_files | tr '\n' '\0' | xargs -0 -n 1 bcftools index --force
    ~{if use_shapeit5_ligate then
      "ligate \\\n" +
      (if cpu > 1 then "  --thread " + cpu + " \\\n" else "") +
      "  --input $vcf_files \\\n" +
      (if defined(pedigree_file) then "  --pedigree \"" + basename(select_first([pedigree_file])) + "\" \\\n" else "") +
      "  --output \"" + filebase + ".bcf\" \\\n" +
      "  --index \\\n" +
      "  1>&2"
    else
      "bcftools concat \\\n" +
      "  --no-version \\\n" +
      "  --output-type b \\\n" +
      "  --output \"" + filebase + ".bcf\" \\\n" +
      "  --compact-PS \\\n" +
      "  --file-list $vcf_files \\\n" +
      "  --ligate \\\n" +
      "  --ligate-force \\\n" +
      "  --write-index"}
    ~{if defined(pedigree_file) then "rm \"" + basename(select_first([pedigree_file])) + "\"" else ""}
    cat $vcf_files | tr '\n' '\0' | xargs -0 rm
    cat $vcf_files | sed 's/$/.csi/' | tr '\n' '\0' | xargs -0 rm
    bcftools query --list-samples "~{filebase}.bcf" | wc -l
  >>>

  output {
    File vcf_file = filebase + ".bcf"
    File vcf_idx = filebase + ".bcf.csi"
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

task vcf_import {
  input {
    File pgt_vcf_file
    File pgt_vcf_idx
    File unphased_vcf_file
    File unphased_vcf_idx

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float pgt_size = size(pgt_vcf_file, "GiB")
  Float vcf_size = size(unphased_vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * pgt_size + 2.0 * vcf_size)])
  Float memory = select_first([memory_override, 3.5])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])
  String filebase = basename(basename(unphased_vcf_file, ".bcf"), ".vcf.gz")

  command <<<
    set -euo pipefail
    mv "~{pgt_vcf_file}" .
    mv "~{pgt_vcf_idx}" .
    mv "~{unphased_vcf_file}" .
    mv "~{unphased_vcf_idx}" .
    bcftools annotate \
      --no-version \
      --output "~{filebase}.phased.bcf" \
      --output-type b \
      --annotations "~{basename(pgt_vcf_file)}" \
      --columns -FMT/GT \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      "~{basename(unphased_vcf_file)}" \
      --write-index
    bcftools query --list-samples "~{filebase}.phased.bcf" | wc -l
    rm "~{basename(pgt_vcf_file)}"
    rm "~{basename(pgt_vcf_idx)}"
    rm "~{basename(unphased_vcf_file)}"
    rm "~{basename(unphased_vcf_idx)}"
  >>>

  output {
    File vcf_file = filebase + ".phased.bcf"
    File vcf_idx = filebase + ".phased.bcf.csi"
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

# uses hack from https://github.com/samtools/bcftools/pull/1505
# the command requires BCFtools 1.14 due to bug https://github.com/samtools/bcftools/issues/1418
task get_max_nrecords {
  input {
    File vcf_idx

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    mv "~{vcf_idx}" .
    bcftools index --stats "~{basename(vcf_idx)}" | cut -f3 | sort -n -T . | tail -n1
    rm "~{basename(vcf_idx)}"
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

# the command requires BCFtools 1.15 or newer due to bug https://github.com/samtools/htslib/issues/1362
task vcf_mocha {
  input {
    Int n_smpls
    Int max_n_markers
    String assembly
    File? rules_file
    File pvcf_file
    File pvcf_idx
    File? sample_tsv_file
    File? xcl_vcf_file
    File? xcl_vcf_idx
    File? cnp_file
    String? mhc_reg
    String? kir_reg
    String? mocha_extra_args
    String ext_string

    String docker
    Int? cpu_override
    Int? disk_size_override
    Float? memory_override
    Int preemptible = 1
    Int maxRetries = 0
  }

  Float pvcf_size = size(pvcf_file, "GiB")
  Float xcl_size = size(xcl_vcf_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + 2.0 * pvcf_size + xcl_size)])
  Float memory = select_first([memory_override, 3.5 + 18.0 * n_smpls * max_n_markers / 1024 / 1024 / 1024])
  Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * ceil(memory / 13) else 1])
  String filebase = basename(basename(basename(basename(pvcf_file, ".bcf"), ".vcf.gz"), '.phased'), "." + ext_string)
  String input_pvcf_file = sub(sub(basename(pvcf_file), "\\." + ext_string + ".bcf$", ".bcf"), "\\." + ext_string + ".vcf.gz$", ".vcf.gz")
  String input_pvcf_idx = sub(sub(sub(basename(pvcf_idx), "\\." + ext_string + "\\.bcf\\.csi$", ".bcf.csi"), "\\." + ext_string + "\\.vcf\\.gz\\.tbi$", ".vcf.gz.tbi"), "\\." + ext_string + "\\.vcf\\.gz\\.csi$", ".vcf.gz.csi")

  command <<<
    set -euo pipefail
    mv "~{pvcf_file}" ~{if basename(pvcf_file) != input_pvcf_file then "\"" + input_pvcf_file + "\"" else "."}
    mv "~{pvcf_idx}" ~{if basename(pvcf_file) != input_pvcf_file then "\"" + input_pvcf_idx + "\"" else "."}
    echo "~{sep("\n", select_all([rules_file, sample_tsv_file, xcl_vcf_file, xcl_vcf_idx, cnp_file]))}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    bcftools +mocha \
      --genome~{if defined(rules_file) then "-file \"" + basename(select_first([rules_file])) else " \"" + assembly}" \
      --no-version \
      ~{if defined(sample_tsv_file) then "--input-stats \"" + basename(select_first([sample_tsv_file])) + "\"" else ""} \
      ~{if defined(xcl_vcf_file) then "--variants \"^" + basename(select_first([xcl_vcf_file])) + "\"" else ""} \
      ~{if defined(cnp_file) then "--cnp \"" + basename(select_first([cnp_file])) + "\"" else ""} \
      ~{if defined(mhc_reg) then "--mhc \"" + mhc_reg + "\"" else ""} \
      ~{if defined(kir_reg) then "--kir \"" + kir_reg + "\"" else ""} \
      ~{if cpu > 1 then "--threads " + (cpu - 1) else ""} \
      --output "~{filebase}.~{ext_string}.bcf" \
      --output-type b \
      --calls "~{filebase}.calls.tsv" \
      --stats "~{filebase}.stats.tsv" \
      --ucsc-bed "~{filebase}.ucsc.bed" \
      "~{input_pvcf_file}" \
      ~{mocha_extra_args} \
      --write-index
    rm "~{input_pvcf_file}"
    rm "~{input_pvcf_idx}"
    echo "~{sep("\n", select_all([rules_file, sample_tsv_file, xcl_vcf_file, xcl_vcf_idx, cnp_file]))}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File mocha_vcf_file = sub(filebase, "\\." + ext_string + "$", "") + "." + ext_string + ".bcf"
    File mocha_vcf_idx = sub(filebase, "\\." + ext_string + "$", "") + "." + ext_string + ".bcf.csi"
    File calls_tsv = filebase + ".calls.tsv"
    File stats_tsv = filebase + ".stats.tsv"
    File ucsc_bed = filebase + ".ucsc.bed"
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

# https://support.terra.bio/hc/en-us/community/posts/360071476431-Terra-fails-to-delocalize-files-listed-through-read-lines-
task mocha_plot {
  input {
    File vcf_file
    File vcf_idx
    File calls_tsv
    File stats_tsv
    File? cyto_file
    Int? n_chrs
    Float? call_rate_thr
    Float baf_auto_thr = 0.03
    String? mocha_plot_extra_args
    Boolean wgs = false
    Boolean do_not_plot_sex_chromosomes = false

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
    echo "~{sep("\n", select_all([vcf_file, vcf_idx, calls_tsv, stats_tsv, cyto_file]))}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    mkdir pngs
    beg_pos=$(head -n1 "~{basename(calls_tsv)}" | tr '\t' '\n' | grep ^beg_)
    end_pos=$(head -n1 "~{basename(calls_tsv)}" | tr '\t' '\n' | grep ^end_)
    awk -F"\t" -v OFS="\t" -v beg_pos=$beg_pos -v end_pos=$end_pos '
      NR==FNR && FNR==1 {for (i=1; i<=NF; i++) f[$i] = i}
      NR==FNR && FNR>1 && (~{if defined(call_rate_thr) then "$(f[\"call_rate\"])<" + select_first([call_rate_thr]) + " || " else ""}$(f["baf_auto"])>~{baf_auto_thr}) {xcl[$(f["sample_id"])]++}
      NR>FNR && FNR==1 {for (i=1; i<=NF; i++) g[$i] = i}
      NR>FNR && FNR>1 {len=$(g["length"]); bdev=$(g["bdev"]); rel_cov=$(g["rel_cov"])}
      NR>FNR && FNR>1 && !($(g["sample_id"]) in xcl) && $(g["type"])!~"^CNP" &&~{if do_not_plot_sex_chromosomes
    then "\n    chrom!=\"X\" && chrom!=\"chrX\" && chrom!=\"Y\" && chrom!=\"chrY\" &&" else ""}
        ( $(g["chrom"])~"X" && $(g["computed_gender"])=="M" || bdev<0.1 || $(g["n50_hets"])<2e5 ) &&
        ( $(g["bdev_se"])!="nan" || $(g["lod_baf_phase"])!="nan" && $(g["lod_baf_phase"]) > 10.0 ) &&
        ( rel_cov<2.1 || bdev<0.05 || len>5e5 && bdev<0.1 && rel_cov<2.5 || len>5e6 && bdev<0.15 ) {
      print $(g["sample_id"]),$(g["chrom"]),$(g[beg_pos]),$(g[end_pos])}' \
      "~{basename(stats_tsv)}" "~{basename(calls_tsv)}" > "~{basename(calls_tsv, ".tsv")}.coords.tsv"
    while IFS=$'\t' read sample_id chrom beg_pos end_pos; do
      mocha_plot.R \
        ~{if defined(cyto_file) then "--cytoband \"" + basename(select_first([cyto_file])) + "\"" else ""} \
        ~{if defined(n_chrs) then "--n-chrs " + select_first([n_chrs]) else ""} \
        ~{if wgs then "--wgs" else ""} \
        --mocha \
        --stats "~{basename(stats_tsv)}" \
        --png "pngs/${sample_id//[:\\\/ $'\t']/_}.${chrom}_${beg_pos}_$end_pos.png" \
        --vcf "~{basename(vcf_file)}" \
        --samples "$sample_id" \
        --regions $chrom:$beg_pos-$end_pos \
        ~{mocha_plot_extra_args}
      echo pngs/${sample_id//[:\\\/ $'\t']/_}.${chrom}_${beg_pos}_$end_pos.png
    done < "~{basename(calls_tsv, ".tsv")}.coords.tsv"
    rm "~{basename(calls_tsv, ".tsv")}.coords.tsv"
    echo "~{sep("\n", select_all([vcf_file, vcf_idx, calls_tsv, stats_tsv, cyto_file]))}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    Directory pngs = "pngs"
    Array[File] png_files = read_lines(stdout())
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

task mocha_summary {
  input {
    File calls_tsv
    File stats_tsv
    Array[File]+ ucsc_beds
    File? cyto_file
    Int? n_chrs
    String filebase
    Float? call_rate_thr
    Float baf_auto_thr = 0.03

    String docker
    Int cpu = 1
    Int disk_size = 10
    Float memory = 3.5
    Int preemptible = 1
    Int maxRetries = 0
  }

  command <<<
    set -euo pipefail
    echo "~{sep("\n", select_all([calls_tsv, stats_tsv, cyto_file]))}" | \
      tr '\n' '\0' | xargs -0 mv -t .
    ucsc_files=~{write_lines(ucsc_beds)}
    ~{if length(ucsc_beds) > 1 then
    "cat $ucsc_files | tr '\\n' '\\0' | xargs -0 mv -t .\n" +
    "sed -i 's/^.*\\///' $ucsc_files\n" +
    "cat $ucsc_files | tr '\\n' '\\0' | xargs -0 cat | \\\n" +
    "  awk '{if ($0~\"^track\") track=$0; else bed[track]=bed[track]$0\"\\n\"}\n" +
    "  END {for (track in bed) printf track\"\\n\"bed[track]}' > \"" + filebase + ".ucsc.bed\"\n" +
    "cat $ucsc_files | tr '\\n' '\\0' | xargs -0 rm"
    else "mv \"" + ucsc_beds[0] + "\" \"" + filebase + ".ucsc.bed\""}
    summary_plot.R \
      ~{if defined(n_chrs) then "--n-chrs " + select_first([n_chrs]) else ""} \
      --stats "~{basename(stats_tsv)}" \
      --calls "~{basename(calls_tsv)}" \
      --call-rate-thr ~{if defined(call_rate_thr) then select_first([call_rate_thr]) else "1.0"} \
      --baf-auto-thr ~{baf_auto_thr} \
      --pdf "~{filebase}.summary.pdf"
    awk -F "\t" 'NR==FNR && FNR==1 {for (i=1; i<=NF; i++) f[$i] = i}
      NR==FNR && FNR>1 && (~{if defined(call_rate_thr) then "$(f[\"call_rate\"])<" + select_first([call_rate_thr]) + " || " else ""}$(f["baf_auto"])>~{baf_auto_thr}) {xcl[$(f["sample_id"])]++}
      NR>FNR && FNR==1 {for (i=1; i<=NF; i++) g[$i] = i; print}
      NR>FNR && FNR>1 {len=$(g["length"]); bdev=$(g["bdev"]); rel_cov=$(g["rel_cov"])}
      NR>FNR && FNR>1 && !($(g["sample_id"]) in xcl) && $(g["type"])!~"^CNP" &&
        ( $(g["chrom"])~"X" && $(g["computed_gender"])=="M" || bdev<0.1 || $(g["n50_hets"])<2e5 ) &&
        ( $(g["bdev_se"])!="nan" || $(g["lod_baf_phase"])!="nan" && $(g["lod_baf_phase"]) > 10.0 ) &&
        ( rel_cov<2.1 || bdev<0.05 || len>5e5 && bdev<0.1 && rel_cov<2.5 || len>5e6 && bdev<0.15 )' \
      "~{basename(stats_tsv)}" "~{basename(calls_tsv)}" > "~{basename(calls_tsv, ".tsv")}.filtered.tsv"
    pileup_plot.R \
      ~{if defined(cyto_file) then "--cytoband \"" + basename(select_first([cyto_file])) + "\"" else ""} \
      ~{if defined(n_chrs) then "--n-chrs " + select_first([n_chrs]) else ""} \
      --stats "~{basename(stats_tsv)}" \
      --calls "~{basename(calls_tsv, ".tsv")}.filtered.tsv" \
      --call-rate-thr ~{if defined(call_rate_thr) then select_first([call_rate_thr]) else "1.0"} \
      --baf-auto-thr ~{baf_auto_thr} \
      --pdf "~{filebase}.pileup.pdf"
    rm "~{basename(calls_tsv, ".tsv")}.filtered.tsv"
    echo "~{sep("\n", select_all([calls_tsv, stats_tsv, cyto_file]))}" | \
      sed 's/^.*\///' | tr '\n' '\0' | xargs -0 rm
  >>>

  output {
    File ucsc_bed = filebase + ".ucsc.bed"
    File summary_pdf = filebase + ".summary.pdf"
    File pileup_pdf = filebase + ".pileup.pdf"
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
