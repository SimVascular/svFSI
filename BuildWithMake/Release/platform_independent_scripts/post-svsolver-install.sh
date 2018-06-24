#/bin/sh
rm -f /usr/local/bin/svsolver-nompi
cp /usr/local/package/REPLACE_SV_VERSION/REPLACE_TIMESTAMP/generic_launch_script /usr/local/bin/svsolver-nompi
chmod a+rx /usr/local/bin/svsolver-nompi
rm -f /usr/local/bin/svsolver-REPLACE_MPI_NAME
cp /usr/local/package/REPLACE_SV_VERSION/REPLACE_TIMESTAMP/generic_launch_script /usr/local/bin/svsolver-REPLACE_MPI_NAME
chmod a+rx /usr/local/bin/svsolver-REPLACE_MPI_NAME
