#!/bin/bash
./assembly.sh \
	&& scp -P 2203 /home/pio/Projekty/workspace/workspace.372/VcfToolbox/target/VcfToolbox-0.0.2-jar-with-dependencies.jar pstawinski@zgm.wum.edu.pl:/mnt/zgmvol/environment/software/VcfToolbox-0.0.2-jar-with-dependencies.jar
