{
inputs = {
 nixpkgs.url = "github:NixOS/nixpkgs/63678e9f3d3afecfeafa0acead6239cdb447574c";
 systems.url = "github:nix-systems/default";
 devenv.url = "github:cachix/devenv";
};

outputs = { self, nixpkgs, devenv, systems, ... } @ inputs:
 let forEachSystem = nixpkgs.lib.genAttrs (import systems); in {
  devShells = forEachSystem (system: let pkgs = nixpkgs.legacyPackages.${system}; in {
   default = devenv.lib.mkShell { inherit inputs pkgs; modules = [{
    packages = with pkgs; [
            git curl google-cloud-sdk duckdb zip
             (python311.withPackages(ps: with ps; [
               pyarrow duckdb h3 rasterio shapely jupyter geopandas matplotlib
               gdal networkx geopy scikit-learn rtree
             ]))
    ];

    scripts.google.exec = ''gsutil -m rsync -avhP gs://open-buildings-data/v3/points_s2_level_4_gzip $1'';
    scripts.parquet.exec = ''time python -W ignore google_to_parquet.py $@'';
    scripts.csv.exec = ''time python -W ignore csv_to_parquet.py $@'';
    scripts.net.exec = ''time python -W ignore netplan/design.py $@'';
    enterShell = ''
     echo "############################################"
     echo "Welcome to the Network Planner"
     echo ""
     echo "google"
     echo "         Download google open buildings point data (42GB)"
     echo "         google ~/data/buildings"
     echo "parquet"
     echo "         Convert to parquet from gzip"
     echo "         parquet  -h3min 5 -h3max 8 -input_dir ~/data/buildings -output_dir ~/data/parquet"
     echo "administrative"
     echo "         Split by country administrative divisions"
     echo "         administrative -country 'UG' -divisions uganda_districts.parquet -column 'district_name' -input_dir ~/data/parquet -output_dir ~/data/country"
     echo "structures"
     echo "         structures ~/data/country/UG/ARUA ~/data/datasets/structures_UG_ARUA.csv"  
     echo "pue"
     echo "         adds a csv file with lat,long, and optional priority(decimal value between 0-1)"
     echo "         pue ~/data/datasets/pue_UG_ARUA.csv ~/data/datasets/structures_UG_ARUA.csv ~/data/datasets/pue_structures_UG_ARUA.csv"   
     echo "net"
     echo "         Runs the two level network design on the points"
     echo "         net -lan_id 604483198646222847 ~/data/h3dist/ ~/data/tlnd/zombo"
     echo "plan"
     echo "         plan ~/data/inputs/UG-ARUA.zip ~/data/outputs/UG-ARUA.zip"
     echo "jupyter"
     echo "         Review the plan with Jupyter notebook"
     echo "         jupyter notebook"
    '';
   }
  ];
 };
 });
};
}
