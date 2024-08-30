#!/bin/bash
NAME=VcfToolbox-0.1.1-jar-with-dependencies.jar

./assembly.sh \
	&& scp -P 2203 ./target/$NAME pstawinski@zgm.wum.edu.pl:/mnt/zgmvol/environment/software/
