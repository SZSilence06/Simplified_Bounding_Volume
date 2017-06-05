#!/bin/bash
GREEN='\033[1;32m'
NC='\033[0m'

print_message()
{
    echo -e "${GREEN}${1}${NC}"
}

print_message "Installing dependencies..."
dependencies=("libcgal-dev" "libboost-all-dev" "libeigen3-dev")
for dep in ${dependencies[@]}
do
    sudo apt-get install ${dep}
done
print_message "Dependencies installed."

build_type=$1
if [[ $build_type != "Debug" ]]
then
    build_type="Release"
fi
print_message "Build type : $build_type."

print_message "Begin building..."
mkdir ${build_type}
cd ${build_type}
cmake -DCMAKE_BUILD_TYPE=${build_type} ../src
make
print_message "Building done."

#copy data and run
print_message "Running..."
cd bin
cp ../../test_data/*.obj .
./sbvgen -s curve1.obj -e 0.04 -r 0.01 -t
   
