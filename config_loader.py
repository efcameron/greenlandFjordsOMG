from pathlib import Path
import yaml

def load_paths():
    base_dir = Path(__file__).resolve().parents[2]
    config_path = Path(__file__).parent / "config.yaml"

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    paths = config["paths"]

    return {
        "base_dir": base_dir,
        "csv_dir": base_dir / paths["csv_dir"],
        "nc_dir": base_dir / paths["nc_dir"],
        "results_dir": base_dir / paths["results_dir"]
    }