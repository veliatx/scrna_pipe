import pandas as pd

from pathlib import Path


def run_scale():
    """
    """
    scale_path = Path('/home/ec2-user/ScaleRna/')
    sample_sheet = Path('inputs/SampleSheet_YL11_scale.csv').absolute()
    run_folder = Path('inputs/230710_VH01198_25_AAAY5CJHV').absolute()
    fastq_dir = Path('inputs/merged_fastq').absolute()

    sample_csv = Path('inputs/samples.csv').absolute()
    genome = Path('inputs/grch38/grch38.json').absolute()
    output_dir = Path('outputs/run_4').absolute()

    nextflow_cmd = (f'nextflow run -bg -profile docker {scale_path} --samples {sample_csv} --samplesCsv {sample_sheet} '
                    f'--genome {genome} --fastqDir {fastq_dir} --outDir {output_dir}')


def run_qc():
    """
    """
    pass