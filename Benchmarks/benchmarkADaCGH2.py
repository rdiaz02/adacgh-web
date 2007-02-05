# -*- coding: iso-8859-15 -*-
"""
Benchmark ADaCGH2 web application.


"""
import time
import unittest
from funkload.FunkLoadTestCase import FunkLoadTestCase
from webunit.utility import Upload


auto_refresh_string = 'This is an autorefreshing page'
MAX_running_time = 3600 * 1 

def common_part_bench(self,                 
                      MAX_running_time = 3600,
                      auto_refresh_string = auto_refresh_string):    
    """ like above, but does not check anything. simply benchmarking"""
    server_url = self.server_url
    
    while True:
        final_body = self.getBody()
        if final_body.find(auto_refresh_string) < 0:
            break
        time.sleep(5)
        self.get(server_url + self.getLastUrl(),
                 description="Get /cgi-bin/checkdone.cgi")
    print 'OK'


class ADaCGH(FunkLoadTestCase):
    """
    This test use a configuration file Adacgh.conf
    (though it ain't really needed).
    """

    def setUp(self):
        """Setting up test."""
        self.logd("setUp")
        self.server_url = 'http://adacgh2.bioinfo.cnio.es'
        ##self.server_url = self.conf_get('main', 'url')


    def testCBS_large(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("large.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'DNAcopy'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        final_output = 'DNA.smooth.region'
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration


    def testHMM_large(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("large.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'HMM'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration


    def testBioHMM_large(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("large.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'BioHMM'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration

    def testCGHseg_large(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("large.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'CGHseg'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration


    def testGLAD_large(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("large.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'GLAD'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration

    def testWavelets_large(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("large.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'Wavelets'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration

    def testACE_large(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("large.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'ACE'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration

    def testPSW_large(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("large.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'PSW'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration






    def testCBS_small(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("small.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'DNAcopy'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        final_output = 'DNA.smooth.region'
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration


    def testHMM_small(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("small.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'HMM'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration


    def testBioHMM_small(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("small.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'BioHMM'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration

    def testCGHseg_small(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("small.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'CGHseg'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration


    def testGLAD_small(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("small.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'GLAD'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration

    def testWavelets_small(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("small.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'Wavelets'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration

    def testACE_small(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("small.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'ACE'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration

    def testPSW_small(self):
        server_url = self.server_url
        self.get(server_url + "/",
            description="Get /")
        self.post(server_url + "/cgi-bin/adacghR.cgi", params=[
            ['acghData', Upload("empty.txt")],
            ['positionInfo', Upload("empty.txt")],
            ['twofiles', 'One.file'],
            ['acghAndPosition', Upload("small.sample")],
            ['centering', 'None'],
            ['methodaCGH', 'PSW'],
            ['Wave.minDiff', '0.25'],
            ['Wave.merge', 'Yes'],
            ['PSW.nIter', '1000'],
            ['PSW.p.crit', '0.15'],
            ['ACE.fdr', '0.15'],
            ['MCR.gapAllowed', '500'],
            ['MCR.alteredLow', '0.03'],
            ['MCR.alteredHigh', '0.97'],
            ['MCR.recurrence', '75'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="CBS; numerical of parallel")
        common_part(self)
        end_time = time.time()
        duration = end_time - start_time
        print duration

    def tearDown(self):
        """Setting up test."""
        self.logd("tearDown.\n")

if __name__ in ('main', '__main__'):
     unittest.main()
