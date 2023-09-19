Cromwell
========

<!--ts-->
   * [Terra](#terra)
   * [Google Cloud Setup](#google-cloud-setup)
   * [Installation](#installation)
      * [mySQL](#mysql)
      * [SLURM](#slurm)
      * [Singularity](#singularity)
      * [Docker](#docker)
   * [Configure Cromwell](#configure-cromwell)
      * [Local Backend](#local-backend)
      * [SLURM Backend](#slurm-backend)
      * [Google Cloud Backend](#google-cloud-backend)
      * [Shared Configuration](#shared-configuration)
   * [Running Cromwell](#running-cromwell)
   * [Debugging Cromwell](#debugging-cromwell)
   * [Troubleshooting](#troubleshooting)
   * [Acknowledgements](#acknowledgements)
<!--te-->

Terra
=====

While WDL pipelines can be run with a Cromwell server alone, even on your laptop, this can be also done in the [Terra](https://terra.bio/) environment developed at the [Broad Institute](https://www.broadinstitute.org/data-sciences-platform) which prevents the user from having to configure the Google Cloud, mySQL and the Cromwell server on his own. To set up you Terra workspace on Terra you will have to:

* [Setup an account](https://support.terra.bio/hc/en-us/articles/360034677651-Account-setup-and-exploring-Terra) with Terra associated with a billing account
* Create a new Terra workspace
* Find the [MoChA method](https://portal.firecloud.org/?return=terra#methods/mocha/mocha/) in the Broad Methods Repository and have it exported to your workspace (you can choose **sample_set** as Root Entity Type)
* Setup resources for GRCh37 or GRCh38 (available for download [here](http://software.broadinstitute.org/software/mocha/)) by unpacking them in a directory in your own Google bucket, such as **gs://{google-bucket}/GRCh38/** and making sure that the **ref_path** variable points to that path
* Once you know the Terra project you will use to run the pipeline, you can use the [SAM API](https://sam.dsde-prod.broadinstitute.org/) identify your Terra service account using the following commands:
```
$ gcloud auth application-default login
$ curl -X GET https://sam.dsde-prod.broadinstitute.org/api/google/v1/user/petServiceAccount/{terra-project} -H Authorization:\ Bearer\ $(gcloud auth application-default print-access-token)
```
* Make sure the bucket with the reference genome resources can be accessed by your Terra service account (see [here](https://support.terra.bio/hc/en-us/articles/360050202631)). Your Terra service account will be labeled like **pet-012345678909876543210@{terra-project}.iam.gserviceaccount.com** and you will have to add it to the list of members with permission to access the bucket with the role "Storage Object Viewer". You can use a command like the following:
```
$ gsutil iam ch serviceAccount:pet-012345678909876543210@{terra-project}.iam.gserviceaccount.com:objectViewer gs://{google-bucket}
```
* Upload your CSV/BPM/EGT/XML/ZIP manifest files in the same location in your Google bucket, such as **gs://{google-bucket}/manifests/** and making sure that the **manifest_path** variable points to that path. Alternatively you can leave the **manifest_path** variable empty and include the full paths in the **batch_tsv_file** table
* Upload your IDAT/GTC/CEL/CHP/TXT/VCF files to your Google bucket, if you have uploaded CHP/TXT files, make sure you upload also the corresponding SNP **AxiomGT1.snp-posteriors.txt** and report **AxiomGT1.report.txt** files
* Format two TSV tables describing samples and batches as explained in the inputs section, upload them to your Google bucket, and make sure variables **sample_tsv_file** and **batch_tsv_file** are set to their location. If you include the file names without their absolute paths, you can include the path column but you have to make sure all data files from the same batch are in the same directory. See below for examples for Illumina and Affymetrix arrays
* Create a JSON file with all the variables including the cohort ID (**sample_set_id**), the **mode** of the analysis, the location of your tables and resource files, parameters to select the amount of desired parallelization. See below for examples for Illumina and Affymetrix arrays
* From your workspace, go to your WORKSPACE tab, select the MoChA workflow that you had previously imported from the Broad Methods Repository
* Select option \"Run workflow with inputs defined by file paths\" and we further recommend to select options \"Use call caching\" and \"Delete intermediate outputs\" (as these can increase 4-5x the storage footprint)
* Click the on \"upload json\" and select the JSON files with the variable defining your run. Then click on \"SAVE\" and then click on \"RUN ANALYSIS\"
* A job will have spawned that you will be able to monitor through the Job Manager to check the pipeline progress. While monitoring the progress, you will be able to open some summary outputs that are generated before the pipeline fully completes
* Notice that the Terra Job Manager has numerous limitations for how you can monitor your workflow. I strongly advise to manually read the metadata of your workflow using the following commands:
```
$ gcloud auth application-default login
$ curl -X GET https://api.firecloud.org/api/workflows/v1/{workflow_uuid}/metadata -H Authorization:\ Bearer\ $(gcloud auth application-default print-access-token) > metadata.{workflow_uuid}.json
```

Google Cloud Setup
==================

If you want to run your own Cromwell server on Google Cloud, you will have to provide you Google service account with the right permissions and set up a virtual machine where to run Cromwell and a mySQL server. These are the steps that I personally advice to take:

* Initialize your Google Cloud configuration:
```
$ gcloud init
```
* If your institution has not provided you with one, create a billing account from the [Billing](https://console.cloud.google.com/billing) page
* Create a project `{google-project}` in your Google Cloud Platform account from the [Manage resources](https://console.cloud.google.com/cloud-resource-manager) page
* Find default service account for `{google-project}` from the [IAM](https://console.cloud.google.com/iam-admin/iam) page which should be labeled as `{google-number}-compute@developer.gserviceaccount.com`
* Download a private json key for your service account through the following command:
```
$ gcloud iam service-accounts keys create {google-project}.key.json --iam-account={google-number}-compute@developer.gserviceaccount.com
```
or from the [Service accounts](https://console.cloud.google.com/iam-admin/serviceaccounts) page by selecting the `{google-project}` first and then the `{google-number}-compute@developer.gserviceaccount.com` service account
* Enable the Cloud Life Sciences API for your Google project from the Google [console](https://console.developers.google.com/apis/library/lifesciences.googleapis.com)
* Enable the Compute Engine API for your Google project from the Google [console](https://console.cloud.google.com/apis/library/compute.googleapis.com)
* Enable the Google Cloud Storage JSON API for your Google project from the Google [console](https://console.cloud.google.com/apis/library/storage-api.googleapis.com)
```
$ for api in lifesciences compute storage-api; do
  gcloud services enable $api.googleapis.com
done
```
* You need the following roles available to your service account: `Cloud Life Sciences Workflows Runner`, `Service Account User`, `Storage Object Admin`. To add these roles you need to have been assigned the rights to change roles you can use the following commands (or manually from the [IAM](https://console.cloud.google.com/iam-admin/iam) page):
```
$ for role in lifesciences.workflowsRunner iam.serviceAccountUser storage.objectAdmin; do
  gcloud projects add-iam-policy-binding {google-project} --member serviceAccount:{google-number}-compute@developer.gserviceaccount.com --role roles/$role
done
```
* Create a Google bucket `gs://{google-bucket}` making sure request pays is not activated from the [Storage](https://console.cloud.google.com/storage/browser) page or through the following command:
```
$ gsutil mb -p {google-project} -l us-central1-a -b on gs://{google-bucket}
```
* Create a Google virtual machine (VM) with name `INSTANCE-ID` from the [VM instances](https://console.cloud.google.com/compute/instances) page. The [n1-standard-1](https://cloud.google.com/compute/vm-instance-pricing#n1_standard_machine_types) (1 vCPU, 3.75 GB memory) VM with Debian GNU/Linux 11 (bullseye) can be sufficient, but make sure to provide at least 200GB of space for the boot disk (switch to Standard persistent disk to limit costs), rather than the default 10GB, as the mySQL database will easily fail with the default space settings. However, to be safe, we advise to select the slightly pricier [N1 custom](https://cloud.google.com/compute/vm-instance-pricing#n1_custommachinetypepricing) with 1 vCPU and 6.5GB memory as Cromwell can require significant amount of memory for large workflows. You can use the following command:
```
$ gcloud compute instances create INSTANCE-ID --project {google-project} --zone us-central1-a --boot-disk-size 200 --boot-disk-type pd-standard --custom-cpu 1 --custom-memory 6.5GB --image-project debian-cloud --image-family debian-11
```

* Start the virtual machine with the following command:
```
$ gcloud compute instances start INSTANCE-ID --project {google-project} --zone us-central1-a
```
* Copy the private json key for your service account to the VM:
```
$ gcloud compute scp {google-project}.key.json INSTANCE-ID: --project {google-project} --zone us-central1-a
```
* Login to the VM with the following command (notice this command will enable ssh tunneling):
```
$ gcloud compute ssh INSTANCE-ID --project {google-project} --zone us-central1-a -- -L 8000:localhost:8000
```
* As the suggested VM costs approximately $35/month to run, when you are not running any computations, to avoid incurring unnecessary costs, you can stop the virtual machine with the following command:
```
$ gcloud compute instances stop INSTANCE-ID --project {google-project} --zone us-central1-a
```

Installation
============

mySQL
-----

* Install the mySQL server on the machine where you plan to run the Cromwell server:
```
$ sudo apt update && sudo apt install default-mysql-server
```
* If you want to configure the directory of the mysql database to something other than `/var/lib/mysql` make sure to set the `datadir` variable in the `/etc/mysql/mysql.conf.d/mysqld.cnf` file
* Start the mySQL server and initialize the root user with the following command (use cromwell as the default root password):
```
$ sudo systemctl start mysql
$ sudo mysql_secure_installation
```
* Login into the mySQL database and run the following commands to create a database to be used by Cromwell (the first two commands might not be needed):
```
$ sudo mysql --user=root --password=cromwell --execute "SET GLOBAL validate_password.policy=LOW"
$ sudo mysql --user=root --password=cromwell --execute "SET GLOBAL validate_password.length=4"
$ sudo mysql --user=root --password=cromwell --execute "CREATE USER 'user'@'localhost' IDENTIFIED BY 'pass'"
$ sudo mysql --user=root --password=cromwell --execute "GRANT ALL PRIVILEGES ON * . * TO 'user'@'localhost'"
$ sudo mysql --user=root --password=cromwell --execute "CREATE DATABASE cromwell"
```

SLURM
-----

* Following the steps [here](https://github.com/broadinstitute/cromwell/blob/develop/src/ci/bin/test_slurm.inc.sh), we can quickly install SLURM on a Debian/Ubuntu based machine:
```
$ sudo apt update && sudo apt install slurm-wlm
```

* Identify available resources on the node that will be running SLURM:
```
$ slurmd -C
```
The output will provide an output such as:
```
NodeName=localhost CPUs=8 Boards=1 SocketsPerBoard=1 CoresPerSocket=4 ThreadsPerCore=2 RealMemory=7419
```
that will be needed for the SLURM configuration file

* To configure SLURM to run on a single node generate the file **slurm.conf** making sure the `Nodename` line matches the output from `slurmd -C` while keeping `State=Unknown` at the end of the line:
```
$ cat << EOF | sudo tee /etc/slurm/slurm.conf > /dev/null
ClusterName=cluster
SlurmctldHost=localhost
MpiDefault=none
ProctrackType=proctrack/linuxproc
ReturnToService=1
SlurmctldPidFile=/run/slurmctld.pid
SlurmdPidFile=/run/slurmd.pid
SlurmdSpoolDir=/var/lib/slurm/slurmd
SlurmUser=slurm
StateSaveLocation=/var/lib/slurm/slurmctld
SwitchType=switch/none
TaskPlugin=task/affinity
SchedulerType=sched/backfill
SelectType=select/cons_tres
SelectTypeParameters=CR_Core_Memory
AccountingStorageType=accounting_storage/none
JobAcctGatherType=jobacct_gather/linux
SlurmctldLogFile=/var/log/slurm/slurmctld.log
SlurmdLogFile=/var/log/slurm/slurmd.log
NodeName=localhost CPUs=8 Boards=1 SocketsPerBoard=1 CoresPerSocket=4 ThreadsPerCore=2 RealMemory=7419 State=UNKNOWN
PartitionName=localpartition Nodes=ALL Default=YES MaxTime=INFINITE State=UP
EOF
```

* To allow CPUs and RAM limits to be correctly managed generate the file **cgroup.conf**:
```
$ cat << EOF | sudo tee /etc/slurm/cgroup.conf > /dev/null
CgroupAutomount=yes
CgroupMountpoint=/sys/fs/cgroup
ConstrainCores=yes
ConstrainDevices=yes
ConstrainKmemSpace=no
ConstrainRAMSpace=yes
ConstrainSwapSpace=yes
EOF
```
Notice that SLURM will not manage the disk space resources necessary to run the tasks (as opposed to when you run Cromwell on Google Cloud) so it is up to you to make sure there is enough space available while running a workflow

* Start the SLURM central management daemon and compute node daemon:
```
$ sudo systemctl start slurmctld
$ sudo systemctl start slurmd
```

Make sure partition and node are in idle state running the following commands:
```
$ sinfo -l && sinfo -Nl
```
You should get an output such as:
```
%a %b %d %H:%M:%S %Y
PARTITION       AVAIL  TIMELIMIT   JOB_SIZE ROOT OVERSUBS     GROUPS  NODES       STATE NODELIST
localpartition*    up   infinite 1-infinite   no       NO        all      1        idle localhost
%a %b %d %H:%M:%S %Y
NODELIST   NODES       PARTITION       STATE CPUS    S:C:T MEMORY TMP_DISK WEIGHT AVAIL_FE REASON
localhost      1 localpartition*        idle 8       1:4:2   7419        0      1   (null) none
```
If instead of **idle** you see **drained**, then you can restore the node with the command:
```
$ sudo scontrol update nodename=localhost state=resume
```

Singularity
-----------

By default Cromwell relies on Dockers but often it is not possible to run Dockers on a shared filesystem. Fortunately it is possible to use [Singularity](https://sylabs.io) as an alternative containerization system. This is provided by Debian/Ubuntu package [singularity-container](https://packages.debian.org/search?keywords=singularity-container) which can be installed either with:

On a Debian machine this can be installed with:
```
$ sudo apt update && sudo apt install singularity-container
```
Or, if the previous does not work, with:
```
$ wget http://ftp.us.debian.org/debian/pool/main/s/singularity-container/singularity-container_3.11.0+ds1-1+b6_amd64.deb
$ sudo apt update && sudo apt install ./singularity-container_3.11.0+ds1-1+b6_amd64.deb
```
Alternatively, if you would rather use [Apptainer](https://apptainer.org), you can use this command instead:
```
$ wget https://github.com/apptainer/apptainer/releases/download/v1.1.8/apptainer_1.1.8_amd64.deb
$ sudo apt update && sudo apt install ./apptainer_1.1.8_amd64.deb
```

When you run Cromwell with an HPC backend, it is possible to download the Dockers with Singularity when they are needed. If you are running Cromwell in a computational environment that does not have internet access, you might have to download the images manually. You can do so with the following command on a node that has internet access:
```
if [ -z $SINGULARITY_CACHEDIR ];
  then CACHE_DIR=$HOME/.singularity/cache
  else CACHE_DIR=$SINGULARITY_CACHEDIR
fi
mkdir -p $CACHE_DIR
DOCKER_REPOSITORY=us.gcr.io/mccarroll-mocha
DOCKER_TAG=1.17-20230919
for docker in debian:stable-slim amancevice/pandas:slim \
  DOCKER_REPOSITORY/{bcftools,autoconvert,iaap_cli,apt,r_mocha,eagle,shapeit4,impute5,beagle5,regenie}:$DOCKER_TAG; do
  DOCKER_NAME=$(sed -e 's/[^A-Za-z0-9._-]/_/g' <<< ${docker})
  IMAGE=$CACHE_DIR/$DOCKER_NAME.sif
  singularity pull $IMAGE docker://${docker}
done
```
You will then have to move manually these images from `$CACHE_DIR` to a location accessible to the nodes that will run the computational tasks. This additional step is only recommended in the rare scenario of no available internet connection available for Cromwell

Docker
------

Docker can be installed with the following command:
```
$ sudo apt update && sudo apt install docker.io
```

It is also possible to install docker without root privileges, though this is only needed if you want to run Cromwell using the **local** backend:
```
$ curl -fsSL https://get.docker.com/rootless | sh
$ systemctl --user start docker
$ export DOCKER_HOST=unix://$XDG_RUNTIME_DIR/docker.sock
```

Configure Cromwell
==================

Once on the node where Cromwell and mySQL will run, install some basic packages as well as the WOMtool and the Cromwell server (replace `XY` with the current version [here](https://github.com/broadinstitute/cromwell/releases)):
```
$ sudo apt update && sudo apt install default-jre-headless jq
$ wget https://github.com/broadinstitute/cromwell/releases/download/XY/womtool-XY.jar
$ wget https://github.com/broadinstitute/cromwell/releases/download/XY/cromwell-XY.jar
```

Local Backend
-------------

If you want to use Cromwell without a job scheduler on a local machine or server, you can use the **local** backend. You will need to create a **cromwell.conf** configuration file including the following backend stanza:
```
backend {
  default = local
  providers {
    local {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        concurrent-job-limit = 4
        run-in-background = true
        filesystems {
            local {
                localization: [
                    "hard-link", "soft-link", "copy"
                ]
                caching {
                    duplication-strategy: [
                        "hard-link", "soft-link", "copy"
                    ]
                    hashing-strategy: "fingerprint"
                    fingerprint-size: 1048576
                    check-sibling-md5: false
                }
            }
        }

        runtime-attributes = """
        String? docker
        String? docker_user
        """

        submit = "/usr/bin/env bash ${script}"

        submit-docker = """
        docker run \
          --rm -i \
          ${"--user " + docker_user} \
          --entrypoint ${job_shell} \
          -v ${cwd}:${docker_cwd} \
          ${docker} ${script}
        """

        root = "cromwell-executions"
        dockerRoot = "/cromwell-executions"
      }
    }
  }
}
```
As Cromwell by default will submit all tasks that are ready to run, it is imperative that you include a [job limit](https://cromwell.readthedocs.io/en/stable/backends/Backends/#backend-job-limits) to make sure the local machine will not run out of memory. Do notice though that by doing so tasks will be executed regardless of whether there is enough memory available. For this reason we advise to use the **local** backend only for testing purposes

SLURM Backend
-------------

Running Cromwell using a High Performance Computing (HPC) framework requires access to a shared filesystem. Because shared filesystem are usually not accessible with root privileges, rather than using Docker images this backend relies on Singularity images so you will need to have Singularity installed. Cromwell can be run on top on several different HPC backends, including [LSF](https://cromwell.readthedocs.io/en/stable/backends/LSF), [SGE](https://cromwell.readthedocs.io/en/stable/backends/SGE), and [SLURM](https://cromwell.readthedocs.io/en/stable/backends/SLURM)

* The main task of setting up a Cromwell server will be to edit the main configuration file. By following the Cromwell documentation (see [here](https://cromwell.readthedocs.io/en/stable/tutorials/Containers/#configuration) and [here](https://cromwell.readthedocs.io/en/stable/tutorials/Containers/#singularity) and [here](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#recognized-runtime-attributes-and-backends) and [here](https://cromwell.readthedocs.io/en/stable/Configuring/#local-filesystem-options)) you will need to create a **cromwell.conf** configuration file including the following backend stanza:
```
backend {
  default = slurm
  providers {
    slurm {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        filesystems {
            local {
                localization: [
                    "hard-link", "soft-link", "copy"
                ]
                caching {
                    duplication-strategy: [
                        "hard-link", "soft-link", "copy"
                    ]
                    hashing-strategy: "fingerprint"
                    fingerprint-size: 1048576
                    check-sibling-md5: false
                }
            }
        }

        runtime-attributes = """
        Int runtime_minutes = 10080
        Int cpu = 1
        Int memory_mb = 3584
        String? docker
        """

        script-epilogue = "sleep 5 && sync"
        submit-docker = """
        # Make sure the SINGULARITY_CACHEDIR variable is set. If not use a default
        # based on the users home.
        if [ -z $SINGULARITY_CACHEDIR ];
            then CACHE_DIR=$HOME/.singularity/cache
            else CACHE_DIR=$SINGULARITY_CACHEDIR
        fi
        # Make sure cache dir exists so lock file can be created by flock
        mkdir -p $CACHE_DIR
        LOCK_FILE=$CACHE_DIR/singularity_pull_flock
        # Create an exclusive filelock with flock. --verbose is useful for
        # debugging, as is the echo command. These show up in `stdout.submit`.
        flock --verbose --exclusive --timeout 900 $LOCK_FILE \
        singularity exec --containall docker://${docker} \
        echo "successfully pulled ${docker}!"

        # Submit the script to SLURM
        sbatch \
          --wait \
          -J=${job_name} \
          -D ${cwd} \
          -o ${out}.sbatch \
          -e ${err}.sbatch \
          -t ${runtime_minutes} \
          -c ${cpu} \
          --mem=${memory_mb} \
          --wrap "singularity exec --containall --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${docker_script}"
        """
        kill-docker = "scancel ${job_id}"
        exit-code-timeout-seconds = 360
        check-alive = "squeue -j ${job_id} || if [[ $? -eq 5004 ]]; then true; else exit $?; fi"
        job-id-regex = "Submitted batch job (\\d+).*"
      }
    }
  }
}
```
If the nodes running the tasks do not have an available internet connection, you might want to remove the section that pulls the Singularity image and replace it with:
```
if [ -z $SINGULARITY_CACHEDIR ];
  then CACHE_DIR=$HOME/.singularity/cache
  else CACHE_DIR=$SINGULARITY_CACHEDIR
fi
DOCKER_NAME=$(sed -e 's/[^A-Za-z0-9._-]/_/g' <<< ${docker})
IMAGE=$CACHE_DIR/$DOCKER_NAME.sif
```
Making sure the Docker images were manually copied in the `$CACHE_DIR` and then replace `docker://${docker}` with `$IMAGE`. The code within the `submit-docker` section will be included in the `script.submit` file that will be located in the `execution` directory of each task

The **filesystems** sub-stanza will instruct Cromwell to activate [fingerprint](https://cromwell.readthedocs.io/en/stable/Configuring/#call-cache-strategy-options-for-local-filesystem) if you are planning to activate CallCaching. This is imperative as on shared filesystems hashes have to be recalculated every time they are requested and this can quickly bring the filesystem down if you are computing hashes of a large number of files. Failure to do so will inevitably lead to **Future Timed Out** or **Futures Timed Out** errors when Cromwell is configured with CallCaching. It is possible that for jobs with many tasks run at once you might still get these type of errors. While the same job can then be restarted using CallCaching, to completely avoid this issue it might be necessary to impose a maximum limit on the [number of tasks](https://cromwell.readthedocs.io/en/stable/backends/Backends/#backend-job-limits) running at the same by setting the `concurrent-job-limit` variable to a reasonable amount. On some NFS setups, when running a job you might get the error `Failed to read_int(".../execution/stdout")` which is the result of asynchronous writing over the `stdout` file between a node running a task and the node running the Cromwell server. To obviate this issue, you might want to include the `script-epilogue = "sleep 5 && sync"` option to give enough time to NFS to synchronozie the `stdout` file across nodes before labeling a task as completed

By setting the [exit-code-timeout-seconds](https://cromwell.readthedocs.io/en/stable/backends/HPC/#exit-code-timeout) option Cromwell will check whether jobs submitted to SLURM are still alive. The `check-alive` command was altered from the default [SLURM](https://cromwell.readthedocs.io/en/stable/backends/SLURM/) configuration file to prevent Cromwell to infer that a task has failed when the command that checks for whether the task is alive times out. The `5004` value corresponds to the SLURM error code `SLURM_PROTOCOL_SOCKET_IMPL_TIMEOUT` as implemented in [slurm_errno.h](https://github.com/SchedMD/slurm/blob/master/slurm/slurm_errno.h) and [slurm_errno.c](https://github.com/SchedMD/slurm/blob/master/src/common/slurm_errno.c)

Google Cloud Backend
--------------------

* Download the [PAPIv2](https://github.com/broadinstitute/cromwell/blob/develop/cromwell.example.backends/PAPIv2.conf) configuration file with the following command:
```
$ wget -O cromwell.conf https://raw.githubusercontent.com/broadinstitute/cromwell/develop/cromwell.example.backends/PAPIv2.conf
```
* This basic configuration file contains only the **backend** configuration stanza. To configure your Cromwell server you have to edit additional configuration stanzas before the **backend** one

* Add the **google** configuration stanza to the `cromwell.conf` configuration file (leave `service-account` and `service_account` verbatim):
```
google {
  application-name = "cromwell"
  auths = [
    {
      name = "service-account"
      scheme = "service_account"
      json-file = "{google-project}.key.json"
    }
  ]
}
```
* Change `auth = "application-default"` to `auth = "service-account"` in the `cromwell.conf` configuration file in both instances where it occurs within the **backend** configuration stanza
* Change `project = "my-cromwell-workflows"` and `project = "google-billing-project"` to `project = "{google-project}"` in the `cromwell.conf` configuration file within the **backend** configuration stanza
* Change `root = "gs://my-cromwell-workflows-bucket"` to `root = "gs://{google-bucket}/cromwell/executions"` in the `cromwell.conf` configuration file within the **backend** configuration stanza
* As Cromwell will need to load some input files to properly organize the batching, it will need the [engine filesystem](https://cromwell.readthedocs.io/en/stable/filesystems/Filesystems/#engine-filesystems) activated for reading files, Add the **engine** configuration stanza to the `cromwell.conf` configuration file (leave `service-account` verbatim):
```
engine {
  filesystems {
    gcs {
      auth = "service-account"
      project = "{google-project}"
    }
  }
}
```

Shared Configuration
--------------------

* Add the **webservice** configuration stanza to the `cromwell.conf` configuration file to make sure the [server](https://cromwell.readthedocs.io/en/stable/Configuring/#server) will only be accessible from the local machine (by default it is open to any interface):
```
webservice {
  port = 8000
  interface = 127.0.0.1
}
```
* Add the following **akka** configuration stanza to the `cromwell.conf` configuration file to make sure Cromwell has enough time to respond to requests for metadata even for large workflows (see [here](https://github.com/hall-lab/cromwell-deployment/blob/master/resources/cromwell-configs/PAPI.v2.conf)), as it would be the case for MoChA runs on biobank-size cohorts:
```
akka {
  http {
    server {
      request-timeout = 300s
      idle-timeout = 300s
    }
  }
}
```
* Add the following **services** configuration stanza to the `cromwell.conf` configuration file to make sure Cromwell can retrieve metadata tables with more than 1,000,000 lines:
```
services {
  MetadataService {
    config {
      metadata-read-row-number-safety-threshold = 10000000
    }
  }
}
```

* Add the **database** configuration stanza to the `cromwell.conf` configuration file (as explained [here](https://cromwell.readthedocs.io/en/stable/Configuring/#database) and [here](https://cromwell.readthedocs.io/en/develop/tutorials/PersistentServer)):
```
database {
  profile = "slick.jdbc.MySQLProfile$"
  db {
    driver = "com.mysql.cj.jdbc.Driver"
    url = "jdbc:mysql://localhost/cromwell?rewriteBatchedStatements=true"
    user = "user"
    password = "pass"
    connectionTimeout = 60000
  }
}
```
* Activate [CallCaching](https://cromwell.readthedocs.io/en/stable/cromwell_features/CallCaching) to allow to re-run jobs without re-running tasks that have already completed by adding the **call-caching** configuration stanza:
```
call-caching {
  enabled = true
  invalidate-bad-cache-results = true
}
```

Running Cromwell
================

* To see if Cromwell can run on your machine, the first step would be to run it on a basic workflow without the use of containers. To do so, edit and generate the following `hello.wdl` workflow file:
```
version development

workflow myWorkflow {
  call myTask
}

task myTask {
  command {
    echo "hello world"
  }
  output {
    String out = read_string(stdout())
  }
}
```
* Run the workflow using Cromwell:
```
$ java -jar cromwell-XY.jar run hello.wdl
```
* You should receive an output as follows:
```
{
  "myWorkflow.myTask.out": "hello world"
}
```
* Once we know that Cromwell can run a basic workflow, start Cromwell as a server on the node running the mySQL server with the following command:
```
$ (java -XX:MaxRAMPercentage=90 -Dconfig.file=cromwell.conf -jar cromwell-XY.jar server &)
```
* If you get the error `Bind failed for TCP channel on endpoint [/127.0.0.1:8000]`, it means some other service is already using port `8000`. Run `lsof -i:8000` to see which service. If it is another Cromwell server, you can stop that with:
```
$ killall java
```
* Edit and generate the following `hello.wdl` workflow file:
```
version development

workflow myWorkflow {
  call myTask
}

task myTask {
  command {
    echo "hello world"
  }
  output {
    String out = read_string(stdout())
  }
  runtime {
    docker: "debian:stable-slim"
  }
}
```
* This workflow will make use of containers so if you are using the shared filesystem backend make sure Cromwell was configured to submit tasks using Singularity. To submit the workflow to the Cromwell server use:
```
$ java -jar cromwell-XY.jar submit hello.wdl
```
* You will receive an output as follows with the `{workflow_uuid}` of the submitted job:
```
[yyyy-mm-dd hh:mm:ss,ss] [info] Slf4jLogger started
[yyyy-mm-dd hh:mm:ss,ss] [info] Workflow 01234567-89ab-dcef-0123-456789abcdef submitted to http://localhost:8000
```
* Create an `options.json` file for Cromwell that should look like this (additional options can be used from [here](https://cromwell.readthedocs.io/en/stable/wf_options/Google)):
```
{
  "delete_intermediate_output_files": true,
  "final_workflow_outputs_dir": "gs://{google-bucket}/cromwell/outputs",
  "use_relative_output_paths": true,
  "final_workflow_log_dir": "gs://{google-bucket}/cromwell/wf_logs",
  "final_call_logs_dir": "gs://{google-bucket}/cromwell/call_logs"
}
```
* Download the MoChA WDL pipeline:
```
$ curl https://raw.githubusercontent.com/freeseek/mocha/master/wdl/mocha.wdl -o mocha.wdl
```
* To verify that your input `{sample-set-id}.json` file is correcly formatted, you can use the WOMtool as follows:
```
$ java -jar womtool-XY.jar validate mocha.wdl -i {sample-set-id}.json
```
* To submit a job, use the command (you can run this straight from your computer if you logged to the VM with ssh tunneling):
```
$ java -jar cromwell-XY.jar submit mocha.wdl -i {sample-set-id}.json -o options.json
```
* If you want to [run](https://cromwell.readthedocs.io/en/stable/CommandLine/#run) a job without starting a server you can use the following syntax instead:
```
$ java -jar cromwell-XY.jar run mocha.wdl -i {sample-set-id}.json -o options.json --metadata-output metadata.json
```
* To read the metadata you can open the file in the Firefox browser:
```
$ firefox metadata.json
```
* To monitor the status of jobs submitted to the server or while running, you can use:
```
$ curl -X GET http://localhost:8000/api/workflows/v1/query | jq ".results[:10][] | {id, name, status, submission}"
```
* To monitor the status of a specific job and extract the metadata:
```
$ curl -X GET http://localhost:8000/api/workflows/v1/{workflow_uuid}/metadata | jq
```
* To monitor whether a job failed due to a missing input file:
```
$ curl -X GET http://localhost:8000/api/workflows/v1/{workflow_uuid}/metadata | jq | grep FileNotFoundException
```
* To extract a summary of the metadata to glance at the workflow progress:
```
$ curl -X GET http://localhost:8000/api/workflows/v1/{workflow_uuid}/metadata | jq '{id, workflowName, labels, status, submission, workflowProcessingEvents, calls: (.calls | map_values(.[-1].executionStatus))}'
```
* To monitor a submitted job with workflow ID `{workflow_uuid}` just open your browser and go to URL:
```
http://localhost:8000/api/workflows/v1/{workflow_uuid}/timing
```
* To abort a running job, you can run:
```
$ curl -X POST http://localhost:8000/api/workflows/v1/{workflow_uuid}/abort | jq
```
* It is possible for the Cromwell server to crash. In this case it is okay to restart the Cromwell server. This will automatically connect with the mySQL server and resume the jobs that were running
* It is also possible for the mySQL database to grow to such an extent that it will then completely fill the VM disk. To check the size of the mySQL database, run:
```
$ sudo ls -l /var/lib/mysql/cromwell
```
* If you want to flush the metadata database to recover disk space from the VM, login into the database and run (though notice that you will have first to stop the Cromwell server and then restart it):
```
$ killall java # to stop the Cromwell server
$ sudo mysql --user=root --password=cromwell --execute "DROP DATABASE cromwell"
$ sudo mysql --user=root --password=cromwell --execute "CREATE DATABASE cromwell"
$ (java -XX:MaxRAMPercentage=90 -Dconfig.file=cromwell.conf -jar cromwell-XY.jar server &)
```
* If you want to clean temporary files and log files, after you have properly moved all the output files you need, you can delete the workflow executions directory defined as `root` in the `cromwell.conf` configuration file. Remember that failure to do proper cleanup could cause you to incur unexpected storage costs. After making sure there are no active jobs running with the Cromwell server, you can remove the executions and logs directory with the following command:
```
$ gsutil -m rm -r gs://{google-bucket}/cromwell/executions gs://{google-bucket}/cromwell/call_logs
```
Notice that either flushing the metadata from the database or removing the workflow executions files will invalidate the cache from previously run tasks

When running computations over large cohorts, it is [best practice](https://cloud.google.com/compute/docs/instances/create-start-preemptible-instance#best_practices) to run the jobs on nights and weekends, as this will increase the chances that preemptible jobs will not be stopped and replaced by a non-preemptible version. This will make jobs run both faster and cheaper

Debugging Cromwell
==================

Whether you run Cromwell on a shared file system or in Google Cloud through Terrra, there are many types of log information that are generated and that can be useful to understand the reasons for a job failing to complete. These can be summarized as: (i) workflow logs; (ii) metadata; and (iii) call logs

While you run a workflow Cromwell will generate a log for the whole submission. This will be output by Cromwell to stderr and, if you have access to the console Cromwell was started from, you will see these logs being generated as the workflow proceeds. These logs can be also located in a file named like `workflow.{workflow_uuid}.log` within the [workflow log directory](https://cromwell.readthedocs.io/en/stable/Configuring/#workflow-log-directory). Issues directly related to the Cromwell server will be located here

To understand what inputs and outputs were provided to each task by the main workflow, the metadata is often the best resource. However, this resource is not automatically generated as a file and you will need to interact with a running Cromwell instance to retrieve it. If you are running Cromwell without mySQL this data will be lost when Cromwell stops running. Assuming Cromwell was configured to respond to the default port `8000`, to retrieve the metadata you will need to run a command such as:
```
$ curl -X GET http://localhost:8000/api/workflows/v1/{workflow_uuid}/metadata > metadata.{workflow_uuid}.json
```
If rather than using Cromwell in server mode you only [run](https://cromwell.readthedocs.io/en/stable/CommandLine/#run) a single workflow, make sure you use the `--metadata-output metadata.json` option to generate the metadata at the end of the run and to visualize the file you can run:
```
$ jq < metadata.{workflow_uuid}.json
```
A common reason for a pipeline failing is that some of the input files were not located where expected. Unfortunately Cromwell does a very poor job at flagging this type of common issues. To see if this is the case, you can grep for the keyword `NoSuchFileException` in the metadata:
```
$ jq < metadata.{workflow_uuid}.json | grep NoSuchFileException
```
For a large cohort the metadata can become very large and it could take some time to extract this information from the Cromwell server. To get a summary of which tasks succeded and failed you can run:
```
$ jq '{id, workflowName, labels, status, submission, workflowProcessingEvents, calls: (.calls | map_values(.[-1].executionStatus))}' < metadata.{workflow_uuid}.json
```

Most of the time when a pipeline fails due to an issue arising within a task, the issue will be best explained within the call logs. If you are running Cromwell using a shared filed system backend, inside the [workflow root directory](https://cromwell.readthedocs.io/en/stable/backends/HPC/#shared-filesystem) you will find files such as:
```
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/script.submit
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/stdout.submit
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/stderr.submit
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/script.check
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/stdout.check
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/stderr.check
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/script
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/stdout
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/stderr
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/rc
```
with `<cromwell_root>` set to `./cromwell-executions` by default. The `rc` file will contain the return code from running the `script` file. If this code is not `0`, the `stderr` files will usually contain important information to figure out what went wrong. The `script.submit` file includes code from the `submit-docker` section of the **cromwell.conf** configuration file. When running Cromwell using a shared filesystem backend it is sometimes possible to attempt to run the commands in the `script.submit` and `script` files to reproduce and understand what caused the errors. The `script.check` file includes code from the `check-alive` section of the **cromwell.conf** configuration file and it is used by Cromwell to check whether the job associated to a given task is still running

If you have the [Run-in-background](https://cromwell.readthedocs.io/en/stable/tutorials/Containers/#cromwell-run-in-background) option `backend.providers.<backend_name>.config.run-in-background = true` in your Cromwell configuration, which is what is recommended if you are running with a shared file system but without a job manager system such as SLURM, then you will find instead:
```
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/script.background
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/stdout.background
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/stderr.background
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/script.submit
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/script
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/stdout
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/stderr
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/execution/rc
```
In this case the `rc` file will contain the return code from the `script.background` file

If you are running Cromwell using a Google Cloud backend, you will find files such as:
```
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/gcs_transfer.sh
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/gcs_localization.sh
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/gcs_delocalization.sh
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/script
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/stdout
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/stderr
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/rc
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/memory_retry_rc
<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>/<call_name>.log
```
The `gcs_transfer.sh` file contains a script of functions required by `gcs_delocalization.sh` and `gcs_localization.sh`. When a VM is spawned to run a task, first `gcs_localization.sh` is run to copy all necessary files within the VM, then docker is run to execute the `script` file, and then the `gcs_delocalization.sh` is run to export the output files to the `<cromwell_root>/<workflow_uuid>/<workflow_root>/call-<call_name>` directory. The output of the docker run will be found in the `stdout` and `stderr` files while the `<call_name>.log` file will include `stderr` and `stdout` of the `gcs_localization.sh` script, the docker run, and the `gcs_delocalization.sh` script. The `memory_retry_rc` file, generated in Cromwell Scala function [checkIfStderrContainsRetryKeys](https://github.com/broadinstitute/cromwell/blob/84/supportedBackends/google/pipelines/common/src/main/scala/cromwell/backend/google/pipelines/common/action/ActionCommands.scala), will contain a value other than `0` if any of the `system.memory-retry-error-keys` strings are observed in the `stderr` file as explained [here](https://cromwell.readthedocs.io/en/stable/cromwell_features/RetryWithMoreMemory/#retry-with-more-memory). Notice that when using a Google Cloud backend there are no `{script,stdour,stderr}.{background,submit}` files as you do not have control of the containerization system to use

While the Cromwell server is running, you can visualize a [timing diagram](https://cromwell.readthedocs.io/en/stable/tutorials/TimingDiagrams/) of the running tasks for a workflow by accessing the URL:
```
http://localhost:8000/api/workflows/v1/{workflow_uuid}/timing
```
If you are running your workflow in Google Cloud, then you will see bars corresponding to the many sub-tasks each task will go through:
* cromwell starting overhead
* Pending
* RequestingExecutionToken
* WaitingForValueStore
* PreparingJob
* CallCacheReading
* RunningJob
* waiting for quota
* Worker "google-pipelines-worker-00112233445566778899aabbccddeeff" assigned in "<runtime_zone>" on a "custom_<cpu>_<memory_mb>" machine
* Pulling "gcr.io/google.com/cloudsdktool/cloud-sdk:354.0.0-alpine"
* Pulling "<docker_image>"
* ContainerSetup
* Background
* Localization
* UserAction
* Delocalization
* Worker released
* Complete in GCE / Cromwell Poll Interval
* UpdatingCallCache
* UpdatingJobStore
* cromwell final overhead

The docker image for each task will be pulled during **Pulling "<docker_image>"**, the `gcs_localization.sh` script will be run during **Localization**, `script` will be run during **UserAction**, and the `gcs_delocalization.sh` script will be run during **Delocalization**

If you cannot figure out on your own what went wrong during a Cromwell run or you receive an error that you do not understand you can contact the [author](mailto:giulio.genovese@gmail.com) but please do your best to first collect and share the workflow logs, the metadata, and all relevant call logs and if the problem seems related to running Cromwell, then include also the Cromwell configuration file you are using

Troubleshooting
===============

These are some of the messages that you might receive when something goes wrong:

* If you run the pipeline on Terra with the `Delete intermediate options` flag selected and your workflow keeps showing as Running even after the final outputs have been generated, it is possible that the Cromwell server behind Terra might have failed while deleting the intermediate outputs. This is an [issue](https://support.terra.bio/hc/en-us/community/posts/360071861791-Job-seems-stuck-indefinitely-at-the-delete-intermediate-files-step-and-does-not-complete) that is being patched
* `Failed to evaluate input 'disk_size' (reason 1 of 1): ValueEvaluator[IdentifierLookup]: No suitable input for ...`: this indicates that Cromwell was unable to find the size of one of the input files for a task, most likely because the file does not exist where indicated by the user
* `The job was stopped before the command finished. PAPI error code 2. Execution failed: generic::unknown: pulling image: docker pull: running ["docker" "pull" "###"]: exit status 1 (standard error: "Error response from daemon: Get https://###: unknown: Project 'project:###' not found or deleted.\n")`: this means that one of the docker images provided does not exist
* `Job exit code 255. Check gs://###/stderr for more information. PAPI error code 9. Please check the log file for more details:`: if this is an error provided by the task cel2chp, it means that the `apt-probeset-genotype` command has encountered an error. Reading the `stderr` file should easily provide an explanation
* `The job was stopped before the command finished. PAPI error code 10. The assigned worker has failed to complete the operation`: this could mean that the job was preempted despite the fact that it was not running in preemptible computing (see [here](https://support.terra.bio/hc/en-us/community/posts/360046714292-Are-you-experiencing-PAPI-error-code-10-Read-this-for-a-workaround-))
* `The compute backend terminated the job. If this termination is unexpected, examine likely causes such as preemption, running out of disk or memory on the compute instance, or exceeding the backend's maximum job duration`: this could be an indication that a task was killed as it requested too much memory
* `The job was aborted from outside Cromwell`: this could be an indication that a task was killed as it requested too much memory
* idat2gtc task fails with internal `stderr` message `Normalization failed! Unable to normalize!`: this means that either the IDAT sample is of very poor quality and it cannot be processed by the GenCall algorithm or that you have matched the IDAT with the wrong Illumina BPM manifest file
* idat2gtc task fails with internall `stderr` message `System.Exception: Unrecoverable Error...Exiting! Unable to find manifest entry ######## in the cluster file!`: this means that you are using the incorrect Illumina EGT cluster file
* If when monitor the status of the job you get the error: `Job Manager is running but encountered a problem getting data from its workflow server. Click here to start over. 500: Internal Server Error` then it means that there is too much metadata input and output into the tasks for the Job Manager to handle te request. This metadata limit is a known issue currently being worked on
* When you mismatch a BPM manifest file with an IDAT in task idat2gtc, iaap-cli, while outputting an error message such as `Normalization failed for sample: ########! This is likely a BPM and IDAT mismatch. ERROR: Index was outside the bounds of the array.` will not return an error code, causing the pipeline to fail at the next task
* If you are running with your own Cromwell server using the PAPIv2 API and you get error `Error attempting to Execute cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager$UserPAPIApiException: Unable to complete PAPI request due to a problem with the request (The caller does not have permission).` whenever Cromwell tries to submit a task, the the cause is the service account that you are using to run the computations with Google Cloud does not have the [Cloud Life Sciences](https://cloud.google.com/life-sciences/docs/concepts/access-control#roles) Workflows Runner (`lifesciences.workflowsRunner`) role set
* If you are running with your own Cromwell server using the PAPIv2 API and some of your tasks start running but they then fail with and you have the error `Error attempting to Execute cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager$UserPAPIApiException: Unable to complete PAPI request due to a problem with the request (Error: checking service account permission: caller does not have access to act as the specified service account: "MY_NUMBER-compute@developer.gserviceaccount.com").` in the workflow log, then the cause is the service account not having the [Service Account User](https://cloud.google.com/iam/docs/service-accounts#user-role) (`roles/iam.serviceAccountUser`) role set
* If you are running with your own Cromwell server using the PAPIv2 API and all of your tasks start running but they all fail with a log file including just the line `yyyy/mm/dd hh:mm:ss Starting container setup.` the cause is the service account not having the [Storage Object](https://cloud.google.com/storage/docs/access-control/iam-roles) Admin (`storage.objectAdmin`) role set
* When running the phasing step on a very large cohort, we have noticed that tasks including the MHC can run very slowly. This is due to the very special nature of the MHC and there is currently no solution. Currently these tasks might dominate the speed of the whole pipeline

Acknowledgements
================

This work is supported by NIH grant [R01 HG006855](http://grantome.com/grant/NIH/R01-HG006855), NIH grant [R01 MH104964](http://grantome.com/grant/NIH/R01-MH104964), NIH grant [R01MH123451](http://grantome.com/grant/NIH/R01-MH123451), and the Stanley Center for Psychiatric Research
