datadir=$1
query=$2

# time calctime 1 &
# time calctime 2 &
# time calctime 3 &
# time calctime 4 &


#./query1 $1 $2 2> /dev/null > q1out &
rm -f endq1 endq1load
./query1 $1 $2 > q1out &
./query2 $1 $2 2> /dev/null > q2out &
./query3 $1 $2 2> /dev/null > q3out &
./query4 $1 $2 > q4out &

#time ./query4 $1 $2 2> /dev/null > q4out &


# ./query1 $1 $2 2> /dev/null > q1out &
# ./query2 $1 $2 2> /dev/null > q2out &
# ./query3 $1 $2 2> /dev/null > q3out &
# ./query4 $1 $2 2> /dev/null > q4out &

#time nice -n 0 ./query4 $1 $2 2>/dev/null > q4out &
#time nice -n 0 ./query4 $1 $2 > q4out &
#time nice -n 0 ./query4 $1 $2 2> q4err > q4out &
#time nice -n 0 ./query1 $1 $2 2> /dev/null > q1out &
#time nice -n 0 ./query1 $1 $2 2> q1err > q1out &
#time ./query2 $1 $2 2> /dev/null > q2out
#time ./query2 $1 $2 2> q2err > q2out
#time ./query3 $1 $2 2> /dev/null > q3out
#time ./query3 $1 $2 2> q3err > q3out

#./query3 $1 $2 2> /dev/null > q3out &

#time ./query4 $1 $2 > q4out &
#./query4 $1 $2 > q4out &
# jobs
wait %1
wait %2
wait %3
touch endq1
wait %4
#sleep 10
rm endq1 endq1load
./main $query q1out q2out q3out q4out
#cat q1err q2err q3err q4err 1>&2
