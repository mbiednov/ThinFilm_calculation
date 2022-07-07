
APPS=pymultilayer

DEPENDENT_LIBS=./multilayer.so

PYTHON_DIR=/home/hongfeng/anaconda3/envs/thinfilm/include/python3.8
NUMPY_DIR=/home/hongfeng/anaconda3/envs/thinfilm/lib/python3.8/site-packages/numpy/core/include/

ALL:${APPS}

%: %.cpp

	swig -python -c++ $@.i
	g++ -c -fPIC -I${PYTHON_DIR} -I${NUMPY_DIR} $^ $@_wrap.cxx
	g++ -shared $@.o $@_wrap.o ${DEPENDENT_LIBS} -O3 -o _$@.so

	rm -r $@.o $@_wrap.o $@_wrap.cxx

	mkdir pymultilayer
	cp _$@.so pymultilayer
	cp $@.py pymultilayer
	cp ${DEPENDENT_LIBS} pymultilayer
	cp test.py pymultilayer

clean:
	rm -r _${APPS}.so
	rm -r $@.o $@_wrap.o $@_wrap.cxx
