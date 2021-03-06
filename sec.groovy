#!/usr/bin/groovy

import org.jlab.io.hipo.HipoDataSource
import org.jlab.io.hipo.HipoDataSync
import org.jlab.detector.base.DetectorType

def reader = new HipoDataSource()
reader.open(args[0])

int n=0
while(n<10) {
  def event = reader.getNextEvent()
  if (event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter")) {
    def evc = event.getBank("REC::Calorimeter")
    def secs = [evc.getShort('pindex')*.toInteger(), evc.getByte('sector')].transpose().collectEntries()
    def evp = event.getBank("REC::Particle")
    evp.getInt("pid").eachWithIndex{pid, ind ->
      println pid+': '+secs[ind]
    }
    n=n+1
  }
}

reader.close()
