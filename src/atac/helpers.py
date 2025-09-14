import pandas as pd
import pyBigWig
import configparser
from pathlib import Path
import os

# --== Config ==--
genome_dir = "datasets/genomes/"
sample_dir = "datasets/Foxk1/"
gtf_file = "gencode.vM36.chr_patch_hapl_scaff.annotation.gff3"
gene_list_file = "datasets/Foxk1/gene_list.txt"
output_file = sample_dir + "target_regions.bed"

bed_file = "datasets/Foxk1/target_regions.bed"
template_ini = "datasets/Foxk1/tfam_study.ini"
output_dir = "datasets/Foxk1/per_region_configs"

# --== GTF/GFF3 Annotation Filtering ==--
def filter_annotations():
  """
  Create a filtered version of the genome annotations to restrict
  which genes are displayed in the final genome track figures.
  """
  
  with open(gene_list_file) as f:
    gene_names = set(line.strip() for line in f if line.strip())
  
  with open(genome_dir + gtf_file) as infile, open(genome_dir + "filtered.gtf", "w") as outfile:
    for line in infile:
      if line.startswith("#"):
        outfile.write(line)
        continue
      # Look for gene_name "XYZ" in attributes column
      if any(f'gene_name "{gene}"' in line for gene in gene_names):
        outfile.write(line)
  
  print(f"Filtered GTF saved to: {output_file}")

# --== Generate Bed File Containing Gene Targets ==--
def generate_bed_file():
  """
  Parse a genome annotation file to collect the regions assocaited
  with each gene of interest and save them to a .bed file
  """
  
  with open(gene_list_file, "r") as f:
    gene_list = set(line.strip() for line in f)
  
  # Read the GFF3 file
  gff3_data = pd.read_csv(genome_dir + gtf_file, sep="\t", header=None, comment="#", 
                          names=["chrom", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"])
  
  genes = gff3_data[gff3_data["feature"] == "gene"]
  
  # Function to extract gene name from the attributes column
  def extract_gene_name(attributes):
      attributes_dict = {}
      for attr in attributes.split(";"):
          key, value = attr.split("=")
          attributes_dict[key] = value
      return attributes_dict.get("gene_name", attributes_dict.get("ID", None))
  
  # Filter genes based on the gene list
  genes["gene_name"] = genes["attributes"].apply(extract_gene_name)
  genes = genes[genes["gene_name"].isin(gene_list)]
  
  # Adjust the start and end positions by 1000 bp
  genes["start"] = genes["start"] - 1000
  genes["end"] = genes["end"] + 1000
  
  # Remove a leading 'chr' from the chrom field
  genes["chrom"] = genes["chrom"].str.replace(r"^chr", "", regex=True)
  
  # Output to BED format
  genes_bed = genes[["chrom", "start", "end", "gene_name"]]
  genes_bed.to_csv(output_file, sep="\t", header=False, index=False)
  
  print(f"Filtered genes saved to {output_file}")

# --== Identify Max Heights & Generate .ini Files ==--
def generate_ini_files(template_ini, output_dir):
  """
  Generate region specific .ini files configured to use
  the maximum peak value across all samples for the figure.
  """
  
  Path(output_dir).mkdir(parents=True, exist_ok=True)

  regions = []
  with open(bed_file) as f:
    for line in f:
      if line.startswith("#") or not line.strip():
        continue
      chrom, start, end, *rest = line.strip().split()
      region_name = rest[0] if rest else f"{chrom}_{start}_{end}"
      regions.append((chrom, int(start), int(end), region_name))
  
  # === Parse .ini template ===
  config = configparser.ConfigParser()
  config.optionxform = str
  config.read(template_ini)
  
  # Identify all bigWig tracks
  bw_tracks = []
  for section in config.sections():
    file_path = config[section].get("file", "")
    if file_path.endswith(".bw") or config[section].get("file_type") == "bigwig":
      bw_tracks.append((section, file_path))

  # === Load bigWig files once ===
  bw_handles = {}
  for section, file_path in bw_tracks:
    if not os.path.isfile(file_path):
      raise FileNotFoundError(f"Missing bigWig file: {file_path}")
    bw_handles[section] = pyBigWig.open(file_path)
  
  # === Process each region ===
  for chrom, start, end, region_name in regions:
      region_config = configparser.ConfigParser()
      region_config.optionxform = str
      region_config.read_dict(config)
  
      # Get the max signal in the region
      region_max = 0
      per_track_max = {}
  
      for section, bw in bw_handles.items():
          try:
              stats = bw.stats(chrom, start, end, type="max")
              val = stats[0] if stats[0] is not None else 0
          except:
              val = 0  # in case region doesn't exist in file
          per_track_max[section] = val
          if val > region_max:
              region_max = val
  
      # Update the config with track-specific max values
      for section in per_track_max:
          region_config[section]["max_value"] = str(region_max)
  
      # Write to new ini
      out_ini = os.path.join(output_dir, f"tracks_{region_name}.ini")
      with open(out_ini, "w") as f:
          region_config.write(f)
  
  print(f"Generated {len(regions)} region-specific .ini files in: {output_dir}")

filter_annotations()
generate_bed_file()
generate_ini_files(sample_dir + "tfam_study.ini",
                   sample_dir + "per_region_configs")
