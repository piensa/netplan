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
     echo "Welcome to the portable STAC catalog example"
     echo ""
     echo "google"
     echo "         Download google open buildings point data (42GB)"
     echo "         google ~/data"
     echo "parquet"
     echo "         Convert to parquet from gzip"
     echo "         parquet ~/data/points_s2_level_4_gzip ~/data/points_parquet/"
     echo "pstac"
     echo "         Create portable stac catalog with data and code"
     echo "         pstac ~/data/google_buildings.zip ~/data/points_parquet"
     echo "jupyter"
     echo "         Edit code and experiment with the data"
     echo "         jupyter notebook"
     echo "git"
     echo "         Version your code"
     echo "         git init ."
     echo "         git add devenv.* *.py"
     echo "         git commit -m 'Initial commit'"
     echo "serve"
     echo "         Serve the data via STAC for apps like QGIS or your colleages on the intranet"
     echo "         serve ~/data/google-buildings.zip"
     echo ""
    '';
   }
  ];
 };
 });
};
}
