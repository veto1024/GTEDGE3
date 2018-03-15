shotnum="118890"
runid="r90"
timeid="1560"
und="_"
for filename in *.txt;
do mv "$filename" "$shotnum$und$timeid$und$runid$und$filename";
done;
