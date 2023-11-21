# netplan
repo for end to end buildings->NP process

# Usage

Get a local copy of the code/scripts.
```bash
git clone https://github.com/SEL-Columbia/netplan
```

Start a shell to run the commands with all the software available.
```
nix --extra-experimental-features 'nix-command flakes' develop --impure
```

Download data from Google Open Buildings (45GB)
```
mkdir ~/data
google ~/data/buildings
```


Run Network
