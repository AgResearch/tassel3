# tassel3

TASSEL is a software package to evaluate traits associations, evolutionary patterns, and linkage disequilibrium. Strengths of this software:

1. It provides a number of new and powerful statistical approaches to association mapping such as a General Linear Model (GLM) and Mixed Linear Model (MLM). MLM is an implementation of the technique which our recently published Nature Genetics paper - [Unified Mixed-Model Method for Association Mapping](https://tassel.bitbucket.io/unified-mixed-model) - which reduces Type I error in association mapping with complex pedigrees, families, founding effects and population structure.

2. Ability to handle a wide range of indels (insertion & deletions). Most software package ignore this type of polymorphism, however, in some species (like maize) this is the most common type of polymorphism.

Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES. (2007) [TASSEL: Software for association mapping of complex traits in diverse samples](https://tassel.bitbucket.io/docs/bradbury2007bioinformatics.pdf). Bioinformatics 23:2633-2635.

TASSEL 3 is no longer being developed by its creators.  This repo contains minor changes for the version in use at AgResearch.

## Documentation

- [TASSEL 3 Genotyping by Sequencing (GBS) pipeline documentation](doc/TasselPipelineGBS.pdf)

## Nix flake

The Nix flake facilitates use of this package either for local development or for incorporation into other Nix flakes.

For example:
```
$ nix develop
$ run_pipeline.pl
Memory Settings: -Xms512m -Xmx1536m
Tassel Pipeline Arguments:
[main] INFO net.maizegenetics.pipeline.TasselPipeline - Tassel Version: 3.0.174  Date: March 3, 2016
[main] WARN net.maizegenetics.pipeline.TasselPipeline - parseArgs: no arguments specified.
```
