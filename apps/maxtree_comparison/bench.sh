#! /bin/bash

BINDIR=/work/carlinet/src/pylene/build/apps/maxtree_comparison


# bench(algo, image, nbits, nthread, size)
# its outputs:
# algo_name, image, nbits, nthread, size, time
function bench
{
    algo=$1
    ima=$2
    nbits=$3
    nthread=$4
    size=$5

    cmdline="$BINDIR/$algo $ima --ntest 10 --sz $size --nthread $nthread --nbits $nbits"
    time=$($cmdline | tee /dev/stderr | grep Run | cut -f 2)
    echo "$algo,$ima,$nbits,$nthread,$size,$time"
}


#SERIALS=(maxtree_serial_berger  maxtree_serial_nister  maxtree_serial_ufindrank
#    maxtree_serial_hqueue  maxtree_serial_ufind_wlc
#    maxtree_serial_najman  maxtree_serial_wilkinson)

SERIALS=(maxtree_serial_najman)



PARALLELS=(maxtree_parallel_hqueue  maxtree_parallel_ufind_line
    maxtree_parallel_pqueue  maxtree_parallel_ufindrank
    maxtree_parallel_ufind)

IMAGES=(/lrde/doc/sip/image/images/samples/barb.pgm
/lrde/doc/sip/image/images/samples/bird.pgm
/lrde/doc/sip/image/images/samples/boat.pgm
/lrde/doc/sip/image/images/samples/bridge.pgm
/lrde/doc/sip/image/images/samples/camera.pgm
/lrde/doc/sip/image/images/samples/frog.pgm
/lrde/doc/sip/image/images/samples/goldhill.pgm
/lrde/doc/sip/image/images/samples/house.pgm
/lrde/doc/sip/image/images/samples/lena.pgm
/lrde/doc/sip/image/images/samples/mandrill.pgm
/lrde/doc/sip/image/images/samples/mountain.pgm
/lrde/doc/sip/image/images/samples/peppers.pgm
/lrde/doc/sip/image/images/samples/washsat.pgm
/lrde/doc/sip/image/images/samples/zelda.pgm
/work/carlinet/img/naturel.pgm
/work/carlinet/img/piscine.pgm
/work/carlinet/img/sat.pgm
/work/carlinet/img/washingtonDC.pgm
)

#SIZES=(1 2 4 8 16 32 64 128 256)
SIZES=(256)
QUANTS=(8 10 12 14 16 18 20 22 24 26 28 30 32)
PQUANTS=(8 10 12 14 16 18)
NTHREADS=(1 2 4 8 12 16)
SIZEREF=8


#first line of the csv
echo "algo_name, serial, image, nbits, nthread, size, time"

########################################
# Compute serial algorithms w.r.t size #
########################################

for algo in ${SERIALS[@]}; do
    for sz in ${SIZES[@]}; do
	for ima in ${IMAGES[@]}; do
	    bench $algo $ima 8 -1 $(($sz * 1024 * 1024))
	done
    done
done
exit 1


########################################
# Compute serial algorithms w.r.t quant #
########################################
for algo in ${SERIALS[@]}; do
    for q in ${QUANTS[@]}; do
	for ima in ${IMAGES[@]}; do
	    bench $algo $ima $q -1 $(($SIZEREF * 1024 * 1024))
	done
    done
done

###########################################
# Compute serial algorithms w.r.t Nthread #
###########################################
for algo in ${PARALLELS[@]}; do
    for th in ${NTHREADS[@]}; do
	for ima in ${IMAGES[@]}; do
	    bench $algo $ima 8 $th $(($SIZEREF * 1024 * 1024))
	done
    done
done

###########################################
# Compute serial algorithms w.r.t Quant #
###########################################
for algo in ${PARALLELS[@]}; do
    for q in ${PQUANTS[@]}; do
	for ima in ${IMAGES[@]}; do
	    bench $algo $ima $q 8 $(($SIZEREF * 1024 * 1024))
	done
    done
done