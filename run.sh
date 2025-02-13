#!/bin/sh

TOP=`pwd`
INFOS_DIRECTORY=$TOP/infos
ROOTFS=$TOP/infos-user/bin/rootfs.tar
KERNEL=$INFOS_DIRECTORY/out/infos-kernel
KERNEL_CMDLINE="boot-device=ata0 init=/usr/init pgalloc.debug=0 pgalloc.algorithm=simple objalloc.debug=0 sched.debug=0 sched.algorithm=cfs syslog=serial $*"
QEMU=/afs/inf.ed.ac.uk/group/teaching/cs3/os/qemu/qemu-3.1.0/x86_64-softmmu/qemu-system-x86_64

$QEMU -kernel $KERNEL -m 4.25G -debugcon stdio -hda $ROOTFS -append "$KERNEL_CMDLINE"
