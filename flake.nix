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
    scripts.lans.exec = ''time python -W ignore csv_to_parquet.py $@'';
    scripts.net.exec = ''time python -W ignore netplan/design.py $@'';
    scripts.plan.exec = ''time python -W ignore plan_network.py $@'';
    enterShell = ''
     echo "############################################"
     echo "Welcome to the Network Planner"
     echo ""
     echo "lans"
     echo "         Convert to parquet from custom CSVs"
     echo "         lans -h3min 6 -h3max 15 -input_dir ~/data/csv -output_dir ~/data/h3dist"
     echo "net"
     echo "         Runs the two level network design on a given lan"
     echo "         net -lan_id 605363480110825471 ~/data/h3dist/ ~/data/tlnd/zombo"
     echo "plan"
     echo "         plan ~/data/tlnd/zombo"
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
