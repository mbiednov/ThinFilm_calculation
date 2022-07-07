APPS = multilayer
LIBAPPS = multilayer.so

all: ${APPS}.o

${APPS}.o: ${APPS}.cpp
	g++ -c $^ -fPIC -o $@
	mv ${APPS}.o python/

${LIBAPPS}: ${APPS}.o
	g++ -O3 -o $@ -shared $^
	rm -r $^
	mv ${LIBAPPS} python/

clean:
	rm -r ${LIBAPPS}
