# mvn clean compile package docker:build && docker push glosbe/wordlist
NAME=genebeam/vcftoolbox

mvn clean package assembly:single


docker build -f Dockerfile -t $NAME . && docker push $NAME
