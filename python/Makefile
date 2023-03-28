
APPS=pymultilayer

DEPENDENT_LIBS=multilayer.o

# COMPILING DEPENDENCY
PYTHON_DIR=/home/hongfeng/anaconda3/envs/thinfilm/include/python3.8/
NUMPY_DIR=/home/hongfeng/anaconda3/envs/thinfilm/lib/python3.8/site-packages/numpy/core/include/

# THE PYTHON LIBRARY YOU WANT TO BE INSTALLED
INSTALL_DIR=/home/hongfeng/anaconda3/envs/thinfilm/lib/python3.8/site-packages/


ALL:${APPS}

%: %.cpp 

	swig -python -c++ $@.i
	g++ -c -fPIC -I${PYTHON_DIR} -I${NUMPY_DIR} $^ $@_wrap.cxx
	g++ -shared $@.o $@_wrap.o ${DEPENDENT_LIBS} -O3 -o _$@.so

	rm -r $@.o $@_wrap.o $@_wrap.cxx

	mkdir pymultilayer
	cp _$@.so pymultilayer
	cp $@.py pymultilayer
#	cp ${DEPENDENT_LIBS} pymultilayer
	cp test.py pymultilayer

	rm -r _${APPS}.so ${APPS}.py multilayer.o

	@echo ${INSTALL_DIR}pymultilayer

	# INSTALLATION
	@if [ -d ${INSTALL_DIR}pymultilayer ]; then\
		rm -r ${INSTALL_DIR}pymultilayer;\
		echo "FOUND: $(INSTALL_DIR)pymultilayer";\
		echo "DELETE: $(INSTALL_DIR)pymultilayer";\
	fi
	cp -r pymultilayer ${INSTALL_DIR}
	rm -r pymultilayer
	# DONE 
	

.PHONY: clean
clean:
	rm -r pymultilayer
	rm -r $@.o $@_wrap.o $@_wrap.cxx
