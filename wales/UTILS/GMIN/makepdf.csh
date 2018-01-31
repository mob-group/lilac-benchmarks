#!/bin/csh
foreach file ([0-9]*.hits.*)
   gminconv < $file > crap ; head -1 crap > pdf.$file
end
rm crap
