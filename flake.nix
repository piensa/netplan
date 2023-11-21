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
    scripts.parquet.exec = ''time python -W ignore google_to_parquet.py $1 $2'';
    scripts.pstac.exec = ''
     zip -r -0 $1 devenv.nix devevn.yaml devenv.lock pstac.zip
     zip -r -0 $1 *.py
     zip -r -0 $1 notebooks/*
     zip -r -0 $1 LICENSE
     cp pstac.zip $1
     zip -r -0 -j $1 $2 
    '';
    scripts.serve.exec = ''
     chmod +x $1
     $1 -m http.serve -d /zip
    '';
    enterShell = ''
     echo "############################################"
     echo "Welcome to the Network Planner"
     echo ""
     echo "google"
     echo "         Download google open buildings point data (42GB)"
     echo "         google ~/data/buildings"
     echo "parquet"
     echo "         Convert to parquet from gzip"
     echo "         parquet -h3min 5 -h3max 8 ~/data/buildings ~/data/parquet_by_country/"
     echo "structures"
     echo "         structures ~/data/parquet_by_country/UG/ARUA ~/data/datasets/structures_UG_ARUA.csv"  
     echo "pue"
     echo "         adds a csv file with lat,long, and optional priority(decimal value between 0-1)"
     echo "         pue ~/data/datasets/pue_UG_ARUA.csv ~/data/datasets/structures_UG_ARUA.csv ~/data/datasets/pue_structures_UG_ARUA.csv"   
     echo "net"
     echo "         Runs the two level network design on the points"
     echo "         net -cutoutmax 2 -cutoutmin 0 ~/data/datasets/pue_structures_UG_ARUA.csv ~/data/datasets/config.json ~/data/inputs/UG-ARUA.zip"
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
