## Get LaPack Repo.
git clone https://github.com/Reference-LAPACK/lapack.git

## Into directory
cd lapack

## Compilation
cp make.inc.example make.inc
make lib -j 2
make blaslib -j 2
