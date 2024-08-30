FROM openjdk:17
ADD target/VcfToolbox-0.1.1-jar-with-dependencies.jar /app.jar

# usage: java -cp $(VCF_TOOLBOX) com.cloudinside.bio.VcfToolbox.JoinVariants ... 
ENTRYPOINT ["/bin/bash"]
