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
            stdenv.mkDerivation rec {
              pname = "tassel3";
              version = "3.0";

              src = fetchzip {
                url = "https://sourceforge.net/projects/tassel/files/Tassel%20${version}/2016-03-03/tassel-${version}-src-20160303.zip";
                hash = "sha256-uoZuvVa4uXXAZkwMUI/W87nH+PBN0kORnPLefonqMOk=";
              };

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

              # unpackPhase = ''
              #   ${unzip}/bin/unzip $src
              # '';

              patches = [
                ./patches/jdk8.patch
                ./patches/unicode.patch
                ./patches/quiet.patch
              ];

              buildPhase = ''
                ant dist
              '';

              installPhase = ''
                runHook preInstall
                mkdir -p $out/bin
                cp run_gbs_pipeline.pl run_pipeline.pl start_tassel.pl $out/bin
                chmod 755 $out/bin/*.pl
                cp -a lib dist $out/bin
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
