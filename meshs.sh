#!/bin/bash

echo "-------Common IDs-------" > MESHs.txt
echo MESH:D007676 >> MESHs.txt 
efetch -db mesh -id D007676 -format native | sed -n '2p' >> MESHs.txt
echo MESH:D007674 >> MESHs.txt 
efetch -db mesh -id D007674 -format native | sed -n '2p' >> MESHs.txt
echo MESH:D007680 >> MESHs.txt 
efetch -db mesh -id D007680 -format native | sed -n '2p' >> MESHs.txt
echo MESH:D007683 >> MESHs.txt 
efetch -db mesh -id D007683 -format native | sed -n '2p' >> MESHs.txt
echo MESH:D058186 >> MESHs.txt 
efetch -db mesh -id D058186 -format native | sed -n '2p' >> MESHs.txt
echo MESH:D007681 >> MESHs.txt 
efetch -db mesh -id D007681 -format native | sed -n '2p' >> MESHs.txt
echo MESH:D007673 >> MESHs.txt 
efetch -db mesh -id D007673 -format native | sed -n '2p' >> MESHs.txt
echo MESH:D007669 >> MESHs.txt 
efetch -db mesh -id D007669 -format native | sed -n '2p' >> MESHs.txt
echo MESH:D052177 >> MESHs.txt 
efetch -db mesh -id D052177 -format native | sed -n '2p' >> MESHs.txt
echo MESH:C537152 >> MESHs.txt 
efetch -db mesh -id C537152 -format native | sed -n '2p' >> MESHs.txt
echo MESH:D007690 >> MESHs.txt 
efetch -db mesh -id D007690 -format native | sed -n '2p' >> MESHs.txt
echo MESH:D021782 >> MESHs.txt 
efetch -db mesh -id D021782 -format native | sed -n '2p' >> MESHs.txt
echo MESH:D012080 >> MESHs.txt 
efetch -db mesh -id D012080 -format native | sed -n '2p' >> MESHs.txt
echo MESH:C538445 >> MESHs.txt 
efetch -db mesh -id C538445 -format native | sed -n '2p' >> MESHs.txt
echo MESH:D000092702 >> MESHs.txt
efetch -db mesh -id D000092702 -format native | sed -n '2p' >> MESHs.txt
echo "-------Unique IDs-------" >> MESHs.txt
echo MESH:D016891 >> MESHs.txt
efetch -db mesh -id D016891 -format native | sed -n '2p' >> MESHs.txt
echo MESH:C531755 >> MESHs.txt
efetch -db mesh -id C531755 -format native | sed -n '2p' >> MESHs.txt
