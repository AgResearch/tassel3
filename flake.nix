{
  description = "Flake for tassel3 that analyzes diversity for sequences, SNPs, or SSRs";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem
      (system:
        let
          pkgs = import nixpkgs {
            inherit system;
          };

          tassel3 =
            with pkgs;
            stdenv.mkDerivation {
              pname = "tassel3";
              version = "3.0";

              src = ./.;

              buildInputs = [
                ant
                jdk8
                perl
              ];

              nativeBuildInputs = [
                perl
                jdk8
                makeWrapper
              ];

              buildPhase = ''
                ant distGBS
              '';

              installPhase = ''
                runHook preInstall
                mkdir -p $out/bin
                cp run_gbs_pipeline.pl run_pipeline.pl start_tassel.pl $out/bin
                chmod 755 $out/bin/*.pl
                cp -a lib $out/bin
                mkdir $out/bin/dist
                cp -a dist/sTASSELGBS.jar $out/bin/dist/sTASSEL.jar
                runHook postInstall
              '';

              postFixup = ''
                for prog in run_gbs_pipeline.pl run_pipeline.pl start_tassel.pl; do
                  wrapProgram $out/bin/''$prog --set PATH "${lib.makeBinPath [
                    jdk8
                  ]}"
                done
              '';
            };

        in
        {
          devShells = {
            default = pkgs.mkShell {
              buildInputs = [ tassel3 ];
            };
          };

          packages = {
            default = tassel3;
          };
        }
      );
}
