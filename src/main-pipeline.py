#!/usr/bin/env python3

import os
import yaml
import logging
from datetime import datetime
from pathlib import Path
import scanpy as sc

from src.preprocessing import preprocess_data, filter_genes
from src.export import export_to_loom
from src.visualization import generate_visualizations
from typing import Dict, Any

def setup_logging(output_dir: str) -> None:
    """Configure logging to both file and console."""
    log_dir = os.path.join(output_dir, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'pipeline_{timestamp}.log')
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def load_config(config_path: str) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config

def create_output_structure(base_dir: str) -> Dict[str, str]:
    """Create output directory structure."""
    directories = {
        'preprocessing': 'preprocessing',
        'pyscenic': 'pyscenic',
        'visualization': 'visualization',
        'logs': 'logs'
    }
    
    paths = {}
    for key, dir_name in directories.items():
        path = os.path.join(base_dir, dir_name)
        os.makedirs(path, exist_ok=True)
        paths[key] = path
    
    return paths

def run_pipeline(config_path: str) -> None:
    """
    Run the complete scRNA-seq analysis pipeline.
    
    Args:
        config_path (str): Path to configuration YAML file
    """
    # Load configuration
    config = load_config(config_path)
    
    # Setup output directory structure
    output_dirs = create_output_structure(config['paths']['output_dir'])
    setup_logging(output_dirs['logs'])
    
    logging.info("Starting scRNA-seq analysis pipeline")
    
    try:
        # 1. Load and preprocess data
        logging.info("Loading data...")
        adata = sc.read(config['paths']['input_data'])
        
        logging.info("Preprocessing data...")
        adata = preprocess_data(
            adata,
            cell_type=config['preprocessing']['cell_type'],
            n_top_genes=config['preprocessing']['n_top_genes']
        )
        
        # Save preprocessed data
        preprocessed_path = os.path.join(
            output_dirs['preprocessing'],
            'preprocessed_data.h5ad'
        )
        adata.write(preprocessed_path)
        
        # 2. Export to loom format
        logging.info("Exporting to loom format...")
        loom_path = os.path.join(output_dirs['pyscenic'], 'output.loom')
        export_to_loom(adata, loom_path)
        
        # 3. Run pySCENIC
        logging.info("Running pySCENIC analysis...")
        scenic_script = os.path.join(os.path.dirname(__file__), 'scripts/run_pyscenic.sh')
        os.system(f"bash {scenic_script} {loom_path} {output_dirs['pyscenic']}")
        
        # 4. Generate visualizations
        logging.info("Generating visualizations...")
        generate_visualizations(
            loom_path=os.path.join(output_dirs['pyscenic'], 'adata_processed_output.loom'),
            adata_path=preprocessed_path,
            output_dir=output_dirs['visualization'],
            group_by=config['visualization']['group_by'],
            top_n_regulons=config['visualization']['top_n_regulons']
        )
        
        logging.info("Pipeline completed successfully!")
        
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}", exc_info=True)
        raise

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Run scRNA-seq analysis pipeline')
    parser.add_argument('--config', required=True, help='Path to config.yaml file')
    args = parser.parse_args()
    
    run_pipeline(args.config)
