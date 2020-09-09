BOOST_LIST="1e0 1e1 1e2 1e3 1e4 1e5"
TRACK_LIST="Quadratic Helix"
SEED_LIST="AB DH"
BOOSTS=$BOOST_LIST SEEDS=$SEED_LIST TRACKS=$TRACK_LIST trackingevalparallel -j 2 -d 5
RESULT=$?
if [ $RESULT -eq 0 ]; then
  exportresults Fittercheck_mfttracks_*.root
  parallel -j 2 boostplots.sh -s={1} -t={2} ::: $SEED_LIST ::: $TRACK_LIST
else
    echo trackingevalparallel Failed!
    return 1
fi

