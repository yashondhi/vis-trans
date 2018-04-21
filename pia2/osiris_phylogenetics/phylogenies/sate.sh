#!/bin/bash
data_file=$1
type=$2
sate_log=$3
auto=$4


if [ "$(auto)" != 'CONFIG' ]; then
/var/www/galaxy_dev/pkgs/satesrc-v2.2.7-2013Feb15/sate-core/run_sate.py -o '.' -d $type -j satejob -i $data_file --auto > $sate_log
else
/var/www/galaxy_dev/pkgs/satesrc-v2.2.7-2013Feb15/sate-core/run_sate.py -o '.' -j satejob -d $type $config_file > $sate_log
fi

cp satejob.marker001.dataset*.dat.aln satejob.aln
