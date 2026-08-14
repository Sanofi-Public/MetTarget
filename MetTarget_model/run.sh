#!/bin/bash

Rscript predict.MetTarget.score.plain.R input.IBD.data.xlsx UC plain.MetTarget.pretrained.UC.RDS > task.UC.log
Rscript predict.MetTarget.score.plain.R input.IBD.data.xlsx CD plain.MetTarget.pretrained.CD.RDS > task.CD.log
Rscript predict.MetTarget.score.plain.R input.IBD.data.xlsx UC.SPHK1 plain.MetTarget.pretrained.UC.RDS > task.UC.SPHK1.log
Rscript predict.MetTarget.score.plain.R input.IBD.data.xlsx CD.SPHK1 plain.MetTarget.pretrained.CD.RDS > task.CD.SPHK1.log
