import csv
import os
import argparse
import requests
from tqdm import tqdm
import math


def fetch_exon_data(gene_name):
    """
    Fetch exon data for a given gene from Ensembl for zebrafish.

    :param gene_name: The name of the gene.
    :return: A list of exons with their start and end positions.
    """
    server = "https://rest.ensembl.org"
    ext = "/lookup/symbol/danio_rerio/{}?expand=1".format(gene_name)

    response = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not response.ok:
        response.raise_for_status()

    data = response.json()
    if 'error' in data:
        return None

    exons = []
    for transcript in data.get('Transcript', []):
        for exon in transcript.get('Exon', []):
            exons.append({'start': exon['start'], 'end': exon['end']})

    return exons


def find_exon(chromosome, coordinate, gene_name):
    """
    Given a chromosome, coordinate, and gene name, return the exon number the coordinate is in.

    :param chromosome: The chromosome where the gene is located.
    :param coordinate: The chromosomal coordinate.
    :param gene_name: The name of the gene.
    :return: The exon number in which the coordinate falls, or None if it doesn't fall within any exon.
    """
    exons = fetch_exon_data(gene_name)
    if not exons:
        return None

    for index, exon in enumerate(exons, start=1):
        if exon['start'] <= coordinate <= exon['end']:
            return index

    return None


def calculate_exon_threshold(exon_count):
    """
    Calculate the exon threshold based on the total number of exons.

    :param exon_count: Total number of exons.
    :return: Exon threshold for filtering.
    """
    if exon_count <= 5:
        # Round up if exon count is 5 or below
        return int(math.ceil(exon_count / 2.0))
    else:
        # Round down if exon count is above 5
        return int(math.floor(exon_count / 2.0))


def read_tsv_files_in_directory_and_write_output(directory_path, output_file_path):
    try:
        # Count the number of TSV files to process for the progress bar
        num_files = sum(1 for file_name in os.listdir(directory_path) if file_name.endswith('.tsv'))

        # Open the output text file for writing
        with open(output_file_path, 'w') as output_file:
            # Initialize the progress bar
            with tqdm(total=num_files, desc="Processing genes") as pbar:
                # Loop through all files in the specified directory
                for file_name in os.listdir(directory_path):
                    if file_name.endswith('.tsv'):
                        tsv_file_path = os.path.join(directory_path, file_name)
                        tsv_file_name = os.path.splitext(file_name)[0]

                        # Extract the gene name from the file name (e.g., "rx3" from "rx3.txt")
                        gene_name = tsv_file_name

                        # Fetch the exons to determine the threshold
                        exons = fetch_exon_data(gene_name)
                        if not exons:
                            print("No exons found for gene {}. Skipping file {}.".format(gene_name, file_name))
                            pbar.update(1)  # Update the progress bar per gene
                            continue

                        # Calculate the exon threshold
                        exon_threshold = calculate_exon_threshold(len(exons))

                        # Open the TSV file and read its content
                        with open(tsv_file_path, 'r') as tsv_file:
                            tsv_reader = csv.reader(tsv_file, delimiter='\t')

                            # Skip the first row
                            next(tsv_reader, None)

                            # Write the header and underline
                            output_file.write(tsv_file_name + '\n')
                            output_file.write('=' * len(tsv_file_name) + '\n')

                            # Read the first 5 valid rows from the TSV file and write them to the text file
                            count = 0
                            for row in tsv_reader:
                                if count >= 5:
                                    break
                                if row:  # Check if row is not empty
                                    # Extract chromosome and position from the third cell
                                    genomic_coordinates = row[2]
                                    chromosome, position = genomic_coordinates.split(':')
                                    chromosome = chromosome.replace('chr', '')
                                    position = int(position)

                                    # Truncate the last three characters of the second cell
                                    second_cell = row[1][:-3] if len(row[1]) > 3 else row[1]

                                    # Find the exon number using the coordinates and gene name
                                    exon_number = find_exon(chromosome, position, gene_name)
                                    if exon_number is not None and exon_number <= exon_threshold:
                                        exon_info = "Exon number: {}".format(exon_number)
                                        # Write to output file
                                        output_file.write(
                                            "[{0}]  Efficiency Score: [{1}]  {2}\n".format(second_cell, row[-1],
                                                                                           exon_info))
                                        count += 1
                            output_file.write("\n")  # Add a newline for separation between files
                        pbar.update(1)  # Update the progress bar per gene
        print("Processed all .tsv files in {} and wrote the output to {}".format(directory_path, output_file_path))
    except Exception as e:
        print("An error occurred: {}".format(e))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Process TSV files in a directory and write the output to a text file.')
    parser.add_argument('-input', type=str, required=True, help='Path to the directory containing .tsv files')
    parser.add_argument('-output', type=str, required=True, help='Path to the output text file')

    args = parser.parse_args()

    read_tsv_files_in_directory_and_write_output(args.input, args.output)