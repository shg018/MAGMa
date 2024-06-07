from django.urls import path
from .views import *

urlpatterns = [
    path('', index, name='index'),
    # path('progressbar', progressbar, name='progressbar'),
    path('home', index, name='home'),
    path('DIALFQ', DIALFQ, name='DIALFQ'),
    path('ProtDIA', ProtDIA, name='ProtDIA'),
    path('PrecDIA', PrecDIA, name='PrecDIA'),
    path('IPDIA', IPDIA, name='IPDIA'),
    path('SILACalltypes', SILACalltypes, name='SILACalltypes'),
    path('SILACIP', SILACIP, name='SILACIP'),
    path('SILACwholeprot', SILACIP, name='SILACwholeprot'),
    path('TMTalltypes', TMTalltypes, name='TMTalltypes'),
    path('TMTfusion', TMTfusion, name='TMTfusion'),
    path('TMT', TMT, name='TMT'),
    path('TMTwholeprot', TMTwholeprot, name='TMTwholeprot'),
    path('TMTphosprot', TMTphosprot, name='TMTphosprot'),
    path('bait_error', bait_error, name='bait_error'),
    path('SN_error', SN_error, name='SN_error'),
    path('intensity_error', intensity_error, name='intensity_error'),
    path('usage', usage, name='usage'),
    path('downloads', downloads, name='downloads'),
    path('about', about, name='about'),
    path('volcanoplot', volcanoplot, name='volcanoplot'),
    path('AdjustedPvalue', AdjustedPvalue, name='AdjustedPvalue'),
    path('FilterForIP', FilterForIP, name='FilterForIP'),
    # path('networkplot', networkplot, name='networkplot'),
    #path(settings.MEDIA_URL, TMT, name='TMT'),
     
]
