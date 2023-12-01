import os
import sys
import glob
import subprocess
import multiprocessing
import pyarrow as pa
import pyarrow.parquet as pq
import duckdb as db
import h3.unstable.vect as vect

def process_file(input):
    filename = input[0]
    params = input[1]
    h3_min = int(params["h3_min"])
    h3_max = int(params["h3_max"])
    out_dir = params["out_dir"]

    print(f"Processing {filename}")
    b = db.sql(f"SELECT latitude, longitude FROM '{filename}'").arrow()
    lats = b["latitude"].combine_chunks().to_numpy()
    lons =  b["longitude"].combine_chunks().to_numpy()
    h3_15 = vect.geo_to_h3(lats, lons, h3_max)
    h3_7 = vect.h3_to_parent(h3_15, h3_min)
    h3_15_arrow = pa.array(h3_15, type=pa.uint64())
    h3_7_arrow = pa.array(h3_7, type=pa.uint64())
    table = pa.table({'h3_max': h3_15_arrow, 'h3_min': h3_7_arrow})
    os.makedirs(out_dir, exist_ok=True)

    target_file = os.path.join(out_dir, os.path.basename(filename).rsplit('.csv.gz', 1)[0]) + ".parquet"
    
    pq.write_table(table, target_file)
    print(f"Wrote {target_file}")

if __name__ == '__main__':
  print(sys.argv)
  # expected format:  -h3min 5 -h3max 8 -input_dir ~/data/buildings -output_dir ~/data/parquet/
  # parsing based on position for now
  params = dict(
    h3_min = sys.argv[2],
    h3_max = sys.argv[4],
    in_dir = sys.argv[6],
    out_dir = sys.argv[8],
  )

  file_paths = glob.glob(f"{params['in_dir']}/*.csv")
  tasks = [(file_path, params) for file_path in file_paths]
  
  print(f"Entering {params['in_dir']}")
  with multiprocessing.Pool(multiprocessing.cpu_count() - 2) as processing_pool:
    processing_pool.map(process_file, tasks)
