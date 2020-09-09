trackingeval -s DH -t Helix -b 1e3
RESULT=$?
if [ $RESULT -eq 0 ]; then
  exportresults Fittercheck_mfttracks_*.root
else
  echo trackingeval failed!
fi
