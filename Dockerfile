FROM java:8
ADD target/VcfToolbox-0.0.2-jar-with-dependencies.jar /app.jar

# usage: java -cp $(VCF_TOOLBOX) com.cloudinside.bio.VcfToolbox.JoinVariants ... 
ENTRYPOINT ["/bin/bash"]
