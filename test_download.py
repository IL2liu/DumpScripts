#!/usr/bin/env python

# DO NOT EDIT THIS FILE BY HAND!
# It is auto-generated from tests.json and gentests.py.

import hashlib
import io
import os
import json
import unittest
import sys

# Allow direct execution
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from youtube_dl.FileDownloader import FileDownloader
import youtube_dl.InfoExtractors

def _file_md5(fn):
    with open(fn, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()
try:
    _skip_unless = unittest.skipUnless
except AttributeError: # Python 2.6
    def _skip_unless(cond, reason='No reason given'):
        def resfunc(f):
            # Start the function name with test to appease nosetests-2.6
            def test_wfunc(*args, **kwargs):
                if cond:
                    return f(*args, **kwargs)
                else:
                    print('Skipped test')
                    return
            return test_wfunc
        return resfunc
_skip = lambda *args, **kwargs: _skip_unless(False, *args, **kwargs)

class DownloadTest(unittest.TestCase):
    PARAMETERS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "parameters.json")

    def setUp(self):
        # Clear old files
        self.tearDown()

        with io.open(self.PARAMETERS_FILE, encoding='utf-8') as pf:
            self.parameters = json.load(pf)

    @_skip_unless(youtube_dl.InfoExtractors.YoutubeIE._WORKING, "IE marked as not _WORKING")
    def test_Youtube(self):
        filename = 'BaW_jenozKc.mp4'
        fd = FileDownloader(self.parameters)
        fd.add_info_extractor(youtube_dl.InfoExtractors.YoutubeIE())
        fd.download(['http://www.youtube.com/watch?v=BaW_jenozKc'])
        self.assertTrue(os.path.exists(filename))
        self.assertEqual(os.path.getsize(filename), 1993883)

    @_skip_unless(youtube_dl.InfoExtractors.DailymotionIE._WORKING, "IE marked as not _WORKING")
    def test_Dailymotion(self):
        filename = 'x33vw9.mp4'
        fd = FileDownloader(self.parameters)
        fd.add_info_extractor(youtube_dl.InfoExtractors.DailymotionIE())
        fd.download(['http://www.dailymotion.com/video/x33vw9_tutoriel-de-youtubeur-dl-des-video_tech'])
        self.assertTrue(os.path.exists(filename))
        md5_for_file = _file_md5(filename)
        self.assertEqual(md5_for_file, 'd363a50e9eb4f22ce90d08d15695bb47')

    @_skip_unless(youtube_dl.InfoExtractors.MetacafeIE._WORKING, "IE marked as not _WORKING")
    def test_Metacafe(self):
        filename = 'aUehQsCQtM.flv'
        fd = FileDownloader(self.parameters)
        fd.add_info_extractor(youtube_dl.InfoExtractors.MetacafeIE())
        fd.add_info_extractor(youtube_dl.InfoExtractors.YoutubeIE())
        fd.download(['http://www.metacafe.com/watch/yt-_aUehQsCQtM/the_electric_company_short_i_pbs_kids_go/'])
        self.assertTrue(os.path.exists(filename))
        self.assertEqual(os.path.getsize(filename), 5754305)

    @_skip_unless(youtube_dl.InfoExtractors.BlipTVIE._WORKING, "IE marked as not _WORKING")
    def test_BlipTV(self):
        filename = '5779306.m4v'
        fd = FileDownloader(self.parameters)
        fd.add_info_extractor(youtube_dl.InfoExtractors.BlipTVIE())
        fd.download(['http://blip.tv/cbr/cbr-exclusive-gotham-city-imposters-bats-vs-jokerz-short-3-5796352'])
        self.assertTrue(os.path.exists(filename))
        md5_for_file = _file_md5(filename)
        self.assertEqual(md5_for_file, '4962f94441605832eb1008eb820ef47a')

    @_skip_unless(youtube_dl.InfoExtractors.XVideosIE._WORKING, "IE marked as not _WORKING")
    def test_XVideos(self):
        filename = '939581.flv'
        fd = FileDownloader(self.parameters)
        fd.add_info_extractor(youtube_dl.InfoExtractors.XVideosIE())
        fd.download(['http://www.xvideos.com/video939581/funny_porns_by_s_-1'])
        self.assertTrue(os.path.exists(filename))
        md5_for_file = _file_md5(filename)
        self.assertEqual(md5_for_file, 'aecab2ea59b7996110a7e409f0c55da3')

    @_skip_unless(youtube_dl.InfoExtractors.VimeoIE._WORKING, "IE marked as not _WORKING")
    @_skip("No output file specified")
    def test_Vimeo(self):
        filename = ''
        fd = FileDownloader(self.parameters)
        fd.add_info_extractor(youtube_dl.InfoExtractors.VimeoIE())
        fd.download(['http://vimeo.com/14160053'])
        self.assertTrue(os.path.exists(filename))
        md5_for_file = _file_md5(filename)
        self.assertEqual(md5_for_file, '1ab4dedc01f771cb2a65e91caa801aaf')

    @_skip_unless(youtube_dl.InfoExtractors.SoundcloudIE._WORKING, "IE marked as not _WORKING")
    def test_Soundcloud(self):
        filename = 'n6FLbx6ZzMiu.mp3'
        fd = FileDownloader(self.parameters)
        fd.add_info_extractor(youtube_dl.InfoExtractors.SoundcloudIE())
        fd.download(['http://soundcloud.com/ethmusic/lostin-powers-she-so-heavy'])
        self.assertTrue(os.path.exists(filename))
        md5_for_file = _file_md5(filename)
        self.assertEqual(md5_for_file, 'c1b9b9ea8bfd620b96b2628664576e1c')

    @_skip_unless(youtube_dl.InfoExtractors.StanfordOpenClassroomIE._WORKING, "IE marked as not _WORKING")
    def test_StanfordOpenClassroom(self):
        filename = 'PracticalUnix_intro-environment.mp4'
        fd = FileDownloader(self.parameters)
        fd.add_info_extractor(youtube_dl.InfoExtractors.StanfordOpenClassroomIE())
        fd.download(['http://openclassroom.stanford.edu/MainFolder/VideoPage.php?course=PracticalUnix&video=intro-environment&speed=100'])
        self.assertTrue(os.path.exists(filename))
        md5_for_file = _file_md5(filename)
        self.assertEqual(md5_for_file, '8aac7873a07dcfaed66b1559ab128514')

    @_skip_unless(youtube_dl.InfoExtractors.CollegeHumorIE._WORKING, "IE marked as not _WORKING")
    @_skip("No output file specified")
    def test_CollegeHumor(self):
        filename = ''
        fd = FileDownloader(self.parameters)
        fd.add_info_extractor(youtube_dl.InfoExtractors.CollegeHumorIE())
        fd.download(['http://www.collegehumor.com/video/6830834/mitt-romney-style-gangnam-style-parody'])
        self.assertTrue(os.path.exists(filename))
        md5_for_file = _file_md5(filename)
        self.assertEqual(md5_for_file, '')

    @_skip_unless(youtube_dl.InfoExtractors.XNXXIE._WORKING, "IE marked as not _WORKING")
    def test_XNXX(self):
        filename = '1135332.flv'
        fd = FileDownloader(self.parameters)
        fd.add_info_extractor(youtube_dl.InfoExtractors.XNXXIE())
        fd.download(['http://video.xnxx.com/video1135332/lida_naked_funny_actress_5_'])
        self.assertTrue(os.path.exists(filename))
        md5_for_file = _file_md5(filename)
        self.assertEqual(md5_for_file, 'c5c67df477eb0d9b058200351448ba4c')

    @_skip_unless(youtube_dl.InfoExtractors.YoukuIE._WORKING, "IE marked as not _WORKING")
    def test_Youku(self):
        filename = 'XNDgyMDQ2NTQw_part00.flv'
        fd = FileDownloader(self.parameters)
        fd.add_info_extractor(youtube_dl.InfoExtractors.YoukuIE())
        fd.download(['http://v.youku.com/v_show/id_XNDgyMDQ2NTQw.html'])
        self.assertTrue(os.path.exists(filename))
        md5_for_file = _file_md5(filename)
        self.assertEqual(md5_for_file, 'ffe3f2e435663dc2d1eea34faeff5b5b')


    def tearDown(self):
        if os.path.exists('BaW_jenozKc.mp4'):
            os.remove('BaW_jenozKc.mp4')
        if os.path.exists('x33vw9.mp4'):
            os.remove('x33vw9.mp4')
        if os.path.exists('aUehQsCQtM.flv'):
            os.remove('aUehQsCQtM.flv')
        if os.path.exists('5779306.m4v'):
            os.remove('5779306.m4v')
        if os.path.exists('939581.flv'):
            os.remove('939581.flv')
        # No file specified for Vimeo
        if os.path.exists('n6FLbx6ZzMiu.mp3'):
            os.remove('n6FLbx6ZzMiu.mp3')
        if os.path.exists('PracticalUnix_intro-environment.mp4'):
            os.remove('PracticalUnix_intro-environment.mp4')
        # No file specified for CollegeHumor
        if os.path.exists('1135332.flv'):
            os.remove('1135332.flv')
        if os.path.exists('XNDgyMDQ2NTQw_part00.flv'):
            os.remove('XNDgyMDQ2NTQw_part00.flv')



if __name__ == '__main__':
    unittest.main()
