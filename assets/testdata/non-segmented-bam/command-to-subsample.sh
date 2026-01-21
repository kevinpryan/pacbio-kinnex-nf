#!/bin/bash -ue
samtools view -b              --subsample 0.0001              --subsample-seed 42              -@ 2              m84039_240124_190648_s3.hifi_reads.bcM0002.bam > m84039_240124_190648_s3_subsampled.bam
