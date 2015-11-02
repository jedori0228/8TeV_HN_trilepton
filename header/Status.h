#include <iostream>
#include <TStopwatch.h>

using namespace std;

void printeventnumber(int index, Long64_t totalevent){
  if(index == 0){
    cout << endl;
  }
  if(index%100000 == 0){
    cout << index << " / " << totalevent << " (" << (double)index*100./(double)totalevent << "%)" << endl;
  }
  if(index == totalevent-1) cout << totalevent << " / " << totalevent << " (100%)" << endl;
}