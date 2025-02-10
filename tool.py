from enum import Enum
from fastapi import FastAPI, HTTPException
from signal import signal, SIGTERM
import sys
import os
import json
import zipfile
from typing import List, Dict, Optional
from pydantic import BaseModel, Field
import tempfile
from pathlib import Path
import logging
import zipfile

# Import core functionality from original implementation
from Bio import SeqIO
import pandas as pd
import requests
import subprocess
from dataclasses import dataclass

# Graceful shutdown for kubernetes
signal(SIGTERM, lambda _1, _2: sys.exit(0))

description = """
A tool for finding orphan genes between sets of species using MMseqs2.
Useful for identifying unique genes that exist in one set of species but not in another.
Input should be two sets of species names and optional parameters.
"""

app = FastAPI(
    title="Unique Gene Finder",
    description=description,
    version=os.environ.get("VERSION", "???"),
    contact={
        "name": "CSIRO",
        "email": "ran12c@csiro.au",
    },
    license_info={
        "name": "MIT",
        "url": "https://opensource.org/license/MIT",
    },
    docs_url="/api",
    root_path=os.environ.get("IVCAP_ROOT_PATH", "")
)

# Core implementation classes from original code


@dataclass
class Species:
    name: str
    tax_id: int
    assembly_accession: Optional[str] = None
    protein_file: Optional[Path] = None


class NCBIDataFetcher:
    """Handle all NCBI API interactions"""

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
        self.headers = {"API_KEY": api_key} if api_key else {}

        # Configure logging
        self.logger = logging.getLogger("NCBIDataFetcher")
        self.logger.handlers = []
        self.logger.propagate = False
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)

        # API endpoints
        self.tax_api = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon_suggest"
        self.assembly_api = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon"
        self.download_api = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"

    def _log_request(self, method: str, url: str, params: Optional[dict] = None, headers: Optional[dict] = None):
        """Log API request details"""
        self.logger.debug(f"\n{'='*80}\nAPI Request:")
        self.logger.debug(f"Method: {method}")
        self.logger.debug(f"URL: {url}")
        if params:
            self.logger.debug(f"Params: {json.dumps(params, indent=2)}")
        if headers:
            safe_headers = {k: v for k,
                            v in headers.items() if k.lower() != 'api_key'}
            self.logger.debug(f"Headers: {json.dumps(safe_headers, indent=2)}")

    def _log_response(self, response: requests.Response):
        """Log API response details"""
        self.logger.debug("\nAPI Response:")
        self.logger.debug(f"Status Code: {response.status_code}")
        self.logger.debug(
            f"Response Time: {response.elapsed.total_seconds():.2f}s")
        try:
            data = response.json()
            json_str = json.dumps(data, indent=2)
            lines = json_str.splitlines()
            if len(lines) > 10:
                truncated = "\n".join(lines[:10]) + \
                    "\n... (response truncated)"
            else:
                truncated = json_str
            self.logger.debug(f"Response Body (truncated):\n{truncated}")
        except:
            self.logger.debug("Response Body: <non-JSON response>")
        self.logger.debug('='*80 + '\n')

    def _make_request(self, method: str, url: str, **kwargs) -> requests.Response:
        """Make HTTP request with logging"""
        self._log_request(method, url, kwargs.get(
            'params'), kwargs.get('headers'))
        response = requests.request(method, url, **kwargs)
        self._log_response(response)
        response.raise_for_status()
        return response

    def get_species_info(self, species_name: str) -> Species:
        """Get taxonomy information for a species"""
        species_name = species_name.replace("_", " ").strip()

        name_map = {
            "escherichia coli": "Escherichia coli",
            "e coli": "Escherichia coli",
            "e. coli": "Escherichia coli",
            "caenorhabditis elegans": "Caenorhabditis elegans",
            "c elegans": "Caenorhabditis elegans",
            "c. elegans": "Caenorhabditis elegans"
        }

        species_name = name_map.get(species_name.lower(), species_name)
        self.logger.info(f"Getting taxonomy info for: {species_name}")

        url = f"{self.tax_api}/{species_name}"
        response = self._make_request('GET', url, headers=self.headers)
        data = response.json()

        if not data or not data.get("sci_name_and_ids"):
            name_parts = species_name.split()
            if len(name_parts) > 2:
                species_name = " ".join(name_parts[:2])
                self.logger.info(
                    f"Retrying with simplified name: {species_name}")
                url = f"{self.tax_api}/{species_name}"
                response = self._make_request('GET', url, headers=self.headers)
                data = response.json()

            if not data or not data.get("sci_name_and_ids"):
                raise ValueError(f"Species not found: {species_name}")

        matches = data["sci_name_and_ids"]
        exact_match = next(
            (m for m in matches if m["sci_name"].lower()
             == species_name.lower()),
            next(
                (m for m in matches if species_name.lower()
                 in m["sci_name"].lower()),
                matches[0]
            )
        )

        tax_id = exact_match["tax_id"]
        official_name = exact_match["sci_name"]
        return Species(name=official_name, tax_id=tax_id)

    def get_assembly_info(self, species: Species) -> Species:
        """Get latest RefSeq assembly for species"""
        self.logger.info(
            f"Getting assembly info for: {species.name} (taxid: {species.tax_id})")
        url = f"{self.assembly_api}/{species.tax_id}/dataset_report"
        params = {
            "filters.assembly_source": "refseq",
            "filters.has_annotation": True,
            "filters.exclude_paired_reports": True,
            "filters.assembly_version": "current"
        }
        response = self._make_request(
            'GET', url, headers=self.headers, params=params)

        data = response.json()
        reports = data.get("reports", [])
        if not reports:
            raise ValueError(f"No assemblies found for {species.name}")

        species.assembly_accession = reports[0]["accession"]
        return species

    def download_proteins(self, species: Species, output_dir: Path) -> Species:
        """Download protein FASTA for species"""
        self.logger.info(
            f"Downloading proteins for: {species.name} (accession: {species.assembly_accession})")
        url = f"{self.download_api}/{species.assembly_accession}/download"
        params = {
            "include_annotation_type": ["PROT_FASTA"],
            "filename": f"{species.name.replace(' ', '_')}_proteins.zip"
        }
        response = self._make_request(
            'GET', url, headers=self.headers, params=params, stream=True)

        zip_file = output_dir / \
            f"{species.name.replace(' ', '_')}_proteins.zip"
        self.logger.debug(f"Saving zip file to: {zip_file}")
        with open(zip_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        output_file = output_dir / \
            f"{species.name.replace(' ', '_')}_proteins.faa"
        self.logger.debug(f"Extracting proteins to: {output_file}")

        try:
            with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                contents = zip_ref.namelist()
                self.logger.debug(f"Zip contents: {contents}")

                protein_files = [
                    name for name in contents if name.endswith('/protein.faa')]
                if not protein_files:
                    raise ValueError(
                        f"No protein.faa file found in zip. Contents: {contents}")

                protein_file = protein_files[0]
                self.logger.debug(f"Found protein file: {protein_file}")

                with zip_ref.open(protein_file) as source, open(output_file, 'wb') as target:
                    target.write(source.read())

                if not output_file.exists():
                    raise ValueError(
                        f"Failed to create output file: {output_file}")
                if output_file.stat().st_size == 0:
                    raise ValueError(f"Extracted file is empty: {output_file}")

                self.logger.debug(
                    f"Successfully extracted {output_file.stat().st_size} bytes")

        except zipfile.BadZipFile:
            raise ValueError(f"Invalid or corrupted zip file: {zip_file}")
        except Exception as e:
            raise ValueError(f"Error extracting protein file: {str(e)}")
        finally:
            if zip_file.exists():
                self.logger.debug(f"Cleaning up zip file: {zip_file}")
                zip_file.unlink()

        species.protein_file = output_file
        return species


class OrphanGeneFinder:
    def __init__(self, set_a_species: List[str], set_b_species: List[str],
                 min_seq_id: float = 0.3, min_coverage: float = 0.8,
                 api_key: Optional[str] = None):
        self.set_a_species = set_a_species
        self.set_b_species = set_b_species
        self.min_seq_id = min_seq_id
        self.min_coverage = min_coverage
        self.ncbi = NCBIDataFetcher(api_key)

        self.logger = logging.getLogger("OrphanGeneFinder")
        self.logger.handlers = []
        self.logger.propagate = False
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)

    def fetch_data(self, temp_dir: Path) -> (List[Species], List[Species]):
        """Fetch all necessary data from NCBI"""
        def process_species_set(species_list: List[str]) -> List[Species]:
            processed = set()
            results = []
            for species_name in species_list:
                species_key = species_name.lower().replace(
                    "_", " ").split("substr.")[0].strip()
                if species_key in processed:
                    continue

                try:
                    species = self.ncbi.get_species_info(species_name)
                    species = self.ncbi.get_assembly_info(species)
                    species = self.ncbi.download_proteins(species, temp_dir)
                    results.append(species)
                    processed.add(species_key)
                except Exception as e:
                    self.logger.error(
                        f"Error processing species {species_name}: {str(e)}")

            return results

        set_a_data = process_species_set(self.set_a_species)
        set_b_data = process_species_set(self.set_b_species)

        if not set_a_data or not set_b_data:
            raise ValueError(
                "Failed to fetch data for one or both species sets")

        return set_a_data, set_b_data

    def _create_mmseqs_db(self, input_fasta: Path, db_path: Path):
        """Create MMseqs2 database"""
        self.logger.info(f"Creating MMseqs2 database for {input_fasta}")
        try:
            result = subprocess.run([
                "mmseqs", "createdb",
                str(input_fasta),
                str(db_path)
            ], check=True, capture_output=True, text=True)
            self.logger.debug(f"MMseqs2 output: {result.stdout}")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"MMseqs2 failed: {e.stderr}")
            raise

    def _run_mmseqs_search(self, query_db: Path, target_db: Path,
                           result_db: Path, tmp_dir: Path):
        """Run MMseqs2 search"""
        self.logger.info("Running MMseqs2 search")
        subprocess.run([
            "mmseqs", "search",
            str(query_db),
            str(target_db),
            str(result_db),
            str(tmp_dir),
            "-s", "7.5",
            "--min-seq-id", str(self.min_seq_id),
            "-c", str(self.min_coverage),
            "--threads", "8",  # Using all available cores
            "--remove-tmp-files", "1"  # Cleanup temp files
        ], check=True, capture_output=True, text=True)

    def _convert_results_to_tsv(self, result_db: Path, query_db: Path, target_db: Path, output_file: Path):
        """Convert MMseqs2 results to TSV format"""
        self.logger.info("Converting results to TSV")
        subprocess.run([
            "mmseqs", "convertalis",
            str(query_db),
            str(target_db),
            str(result_db),
            str(output_file),
            "--format-output", "query,target,pident,qcov,tcov,evalue"
        ], check=True, capture_output=True, text=True)

    def _identify_orphans(self, hits_file: Path, query_genes: Dict) -> pd.DataFrame:
        """Identify orphan genes from MMseqs2 results"""
        self.logger.info("Identifying orphan genes")

        if hits_file.exists() and hits_file.stat().st_size > 0:
            hits = pd.read_csv(hits_file, sep='\t',
                               names=['query', 'target', 'pident', 'qcov', 'tcov', 'evalue'])
            genes_with_hits = set(hits['query'])
        else:
            genes_with_hits = set()

        orphans = [
            {
                'Gene_ID': gene_id.split('|')[1],
                'Species': gene_id.split('|')[0].replace('_', ' '),
                'Length_AA': data['length'],
                'Description': data['description']
            }
            for gene_id, data in query_genes.items()
            if gene_id not in genes_with_hits
        ]

        return pd.DataFrame(orphans)

    def generate_markdown_report(self, orphans_df: pd.DataFrame, total_genes: int) -> str:
        """Generate markdown report of orphan genes"""
        orphan_count = len(orphans_df)
        percentage = (orphan_count/total_genes)*100 if total_genes > 0 else 0

        top_50_orphans = orphans_df.nlargest(50, 'Length_AA')

        report = [
            f"# Orphan Genes Analysis Report\n",
            "## Analysis Parameters",
            f"- Query Species: {', '.join(self.set_a_species)}",
            f"- Target Species: {', '.join(self.set_b_species)}",
            f"- Minimum Sequence Identity: {self.min_seq_id}",
            f"- Minimum Coverage: {self.min_coverage}\n",
            "## Results Summary",
            f"- Total genes analyzed: {total_genes:,}",
            f"- Total orphan genes found: {orphan_count:,}",
            f"- Percentage orphan: {percentage:.1f}%\n",
            "## Top 50 Longest Orphan Genes",
            top_50_orphans.to_markdown(
                index=False) if not orphans_df.empty else "No orphan genes found."
        ]

        return "\n".join(report)

    def find_orphans(self) -> Dict:
        """Find orphan genes using MMseqs2"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            work_dir = Path(tmp_dir)
            set_a_data, set_b_data = self.fetch_data(work_dir)

            if not set_a_data:
                raise ValueError("No Set A data available for analysis")
            if not set_b_data:
                raise ValueError("No Set B data available for analysis")

            query_genes = {}
            with open(work_dir / "combined_query.faa", 'w') as out_f:
                for sp in set_a_data:
                    for record in SeqIO.parse(sp.protein_file, "fasta"):
                        record.id = f"{sp.name.replace(' ', '_')}|{record.id}"
                        record.description = ''
                        query_genes[record.id] = {
                            'length': len(record.seq),
                            'description': ''
                        }
                        SeqIO.write(record, out_f, "fasta")

            total_genes = len(query_genes)

            with open(work_dir / "combined_targets.faa", 'w') as out_f:
                for species in set_b_data:
                    for record in SeqIO.parse(species.protein_file, "fasta"):
                        record.id = f"{species.name.replace(' ', '_')}|{record.id}"
                        record.description = ''
                        SeqIO.write(record, out_f, "fasta")

            query_db = work_dir / "query_db"
            target_db = work_dir / "target_db"
            result_db = work_dir / "result_db"
            hits_file = work_dir / "hits.tsv"

            self._create_mmseqs_db(work_dir / "combined_query.faa", query_db)
            self._create_mmseqs_db(
                work_dir / "combined_targets.faa", target_db)
            self._run_mmseqs_search(query_db, target_db, result_db, work_dir)
            self._convert_results_to_tsv(
                result_db, query_db, target_db, hits_file)

            orphans_df = self._identify_orphans(hits_file, query_genes)
            report_md = self.generate_markdown_report(orphans_df, total_genes)

            top_50_orphans = orphans_df.nlargest(
                50, 'Length_AA').to_dict(orient="records")

            return {
                "total_genes": total_genes,
                "orphan_count": len(orphans_df),
                "percentage_orphan": (len(orphans_df) / total_genes * 100.0) if total_genes > 0 else 0.0,
                "orphans": top_50_orphans,
                "report_markdown": report_md
            }

# FastAPI service models


class ServiceProps(BaseModel):
    min_seq_id: float = Field(
        description="Minimum sequence identity (0.0-1.0)", default=0.3)
    min_coverage: float = Field(
        description="Minimum coverage (0.0-1.0)", default=0.8)
    api_key: Optional[str] = Field(
        description="Optional NCBI API key", default=None)


class ActionProps(BaseModel):
    set_a_species: List[str] = Field(
        description="List of species names to find orphans in")
    set_b_species: List[str] = Field(
        description="List of species names to compare against")


class Props(BaseModel):
    action: ActionProps
    service: ServiceProps


class Response(BaseModel):
    result: str
    report: str


@app.get("/")
def info():
    return {
        "$schema": "urn:sd.platform:schema:ai-tool.1",
        "name": "unique_gene_finder",
        "description": description,
        "action_schema": ActionProps.model_json_schema(by_alias=False),
        "service_schema": ServiceProps.model_json_schema(),
    }


@app.post("/")
def find_orphans(req: Props) -> Response:
    """Find orphan genes between two sets of species"""
    try:
        finder = OrphanGeneFinder(
            set_a_species=req.action.set_a_species,
            set_b_species=req.action.set_b_species,
            min_seq_id=req.service.min_seq_id,
            min_coverage=req.service.min_coverage,
            api_key=req.service.api_key
        )

        result = finder.find_orphans()

        return Response(
            result=json.dumps(result["orphans"]),
            report=result["report_markdown"]
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/_healtz")
def healtz():
    return {"version": os.environ.get("VERSION", "???")}
