{
  description = "cluster-coefficients";

  inputs.flake-utils.url = "github:numtide/flake-utils";
  inputs.nixpkgs.url = "github:nixos/nixpkgs/nixos-24.05";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem
      (system: let
        pkgs = nixpkgs.legacyPackages.${system};
      in
        rec {
          packages = {
            cluster-coefficients = with pkgs; stdenv.mkDerivation {
              name = "cluster-coefficients";
              version = "0.1.0";
              src = lib.sources.cleanSource ./.;

              nativeBuildInputs = [ autoreconfHook pkg-config ];
              buildInputs = [ ginac boost ];
            };
          };

          defaultPackage = packages.cluster-coefficients;

          devShell = with pkgs; mkShell {
            nativeBuildInputs = [ pkg-config autoconf automake ];
            buildInputs = [ ginac boost ];
          };
        }
      );
}
