APPS = multilayer
LIBAPPS = multilayer.so

all: ${LIBAPPS}

${APPS}.o: ${APPS}.cpp
	g++ -c $^ -fPIC -o $@

${LIBAPPS}: ${APPS}.o
	g++ -O3 -o $@ -shared $^
	rm -r $^
	mv ${LIBAPPS} python/

clean:
	rm -r ${LIBAPPS}
