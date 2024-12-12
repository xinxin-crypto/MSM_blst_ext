all: main_test

# One need to modify the following two directories to the correct paths.
openssl_include_directory = /opt/homebrew/include/
openssl_library_directory = /opt/homebrew/Cellar/openssl@3/3.2.1/lib

ifeq ($(group), 1) 
    objects = main_p1.cpp
	output = main_test_p1
else ifeq ($(group), 2)
    objects = main_p2.cpp
	output = main_test_p2
endif

selected_config = ches_config_files/config_file_n_exp_$(config).h

set_config:
	cp $(selected_config) ches_config_files/config_file.h

libblst.a:
	./build.sh

main_test: $(objects) libblst.a set_config
	g++ -std=c++17 -o $(output) -g -O2 $(objects) libblst.a -I$(openssl_include_directory) -L$(openssl_library_directory) -lcrypto

clean:
	rm -f libblst.a
	rm -f main_test_p1
	rm -f main_test_p2