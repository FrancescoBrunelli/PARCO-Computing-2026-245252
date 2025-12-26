#!/bin/bash
wget https://suitesparse-collection-website.herokuapp.com/MM/Hohn/fd12.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Hohn/sinc12.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Averous/epb1.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Schenk_IBMNA/c-46.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Belcastro/human_gene2.tar.gz

tar -xzf fd12.tar.gz
tar -xzf sinc12.tar.gz
tar -xzf c-46.tar.gz
tar -xzf epb1.tar.gz
tar -xzf human_gene2.tar.gz

mv fd12/fd12.mtx .
mv sinc12/sinc12.mtx .
mv c-46/c-46.mtx .
mv epb1/epb1.mtx .
mv human_gene2/human_gene2.mtx .

rm -rf fd12.tar.gz fd12/ sinc12.tar.gz sinc12/ c-46.tar.gz c-46/ epb1.tar.gz epb1/ human_gene2.tar.gz human_gene2/
