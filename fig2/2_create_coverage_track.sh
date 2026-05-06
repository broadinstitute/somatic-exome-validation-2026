#!/bin/bash

python3 scripts/SparK.py \
	-pr 19:36208481-36211352 \
	-o chr19_exp1 \
	-gtf gencode.v19.annotation.gtf \
	-gl WES \
	-l V6 V2 \
	-cg 1 1 1 \
	-tg 1 1 1 \
	-w 100 \
	-tf ICE_GBM.ICB-131-16.577.2-1.u.Pre.bedgraph ICE_GBM.ICB-36-16.111.B3.Pre.bedgraph ICE_GBM.ICB-7-13.971.A1.Pre.bedgraph \
        -cf TWIST_GBM.ICB-131-16.577.2-1.u.Pre.bedgraph TWIST_GBM.ICB-36-16.111.B3.Pre.bedgraph TWIST_GBM.ICB-7-13.971.A1.Pre.bedgraph \
	-pt STD

python3 scripts/SparK.py \
	-pr 1:27022197-27024597 \
	-o chr1_exp3 \
	-gtf gencode.v19.annotation.gtf \
	-gl WES \
	-l V6 V2 \
	-cg 1 1 1 \
	-tg 1 1 1 \
	-w 100 \
	-tf ICE_GBM.ICB-131-16.577.2-1.u.Pre.bedgraph ICE_GBM.ICB-36-16.111.B3.Pre.bedgraph ICE_GBM.ICB-7-13.971.A1.Pre.bedgraph \
        -cf TWIST_GBM.ICB-131-16.577.2-1.u.Pre.bedgraph TWIST_GBM.ICB-36-16.111.B3.Pre.bedgraph TWIST_GBM.ICB-7-13.971.A1.Pre.bedgraph \
	-pt STD

python3 scripts/SparK.py \
	-pr 7:116338189-116341126 \
	-o chr7_exp2 \
	-gtf gencode.v19.annotation.gtf \
	-gl WES \
	-l V6 V2 \
	-cg 1 1 1 \
	-tg 1 1 1 \
	-w 100 \
	-tf ICE_GBM.ICB-131-16.577.2-1.u.Pre.bedgraph ICE_GBM.ICB-36-16.111.B3.Pre.bedgraph ICE_GBM.ICB-7-13.971.A1.Pre.bedgraph \
        -cf TWIST_GBM.ICB-131-16.577.2-1.u.Pre.bedgraph TWIST_GBM.ICB-36-16.111.B3.Pre.bedgraph TWIST_GBM.ICB-7-13.971.A1.Pre.bedgraph \
	-pt STD
