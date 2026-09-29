from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
pdfmetrics.registerFont(TTFont("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf"))
pdfmetrics.registerFont(TTFont("Arial-Bold", "/System/Library/Fonts/Supplemental/Arial Bold.ttf"))
from pathlib import Path
import html
out=Path('output/pdf'); W,H=1000,1210
scale=180/25.4*72/W
c=canvas.Canvas(str(out/'dengue_study_workflow.pdf'),pagesize=(W*scale,H*scale),initialFontName="Arial")
c.setTitle('Workflow for Bayesian modeling of municipal dengue incidence')
c.setAuthor('')
c.scale(scale,scale)
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="217.8mm" viewBox="0 0 {W} {H}">','<rect width="1000" height="1210" fill="white"/>']
ink='#252A30'; colors={'data':('#EDF3F8','#56758A'),'prep':('#FFFFFF','#737B82'),'model':('#FFFFFF','#737B82'),'eval':('#EDF5F2','#527B70')}
def text(x,y,t,size=17,bold=False,anchor='middle',color=ink):
    font='Arial-Bold' if bold else 'Arial'
    c.setFont(font,size); c.setFillColor(HexColor(color))
    {'middle':c.drawCentredString,'start':c.drawString,'end':c.drawRightString}[anchor](x,H-y,t)
    svg.append(f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" font-weight="{700 if bold else 400}" text-anchor="{anchor}" fill="{color}">{html.escape(t)}</text>')
def line(points,color=ink,width=1.5,arrow=False):
    c.setStrokeColor(HexColor(color)); c.setLineWidth(width)
    p=c.beginPath();p.moveTo(points[0][0],H-points[0][1])
    for x,y in points[1:]:p.lineTo(x,H-y)
    c.drawPath(p)
    svg.append(f'<polyline points="{" ".join(f"{x},{y}" for x,y in points)}" fill="none" stroke="{color}" stroke-width="{width}"/>')
    if arrow:
        import math
        x,y=points[-1];px,py=points[-2];a=math.atan2(y-py,x-px)
        pts=[(x,y),(x-8*math.cos(a)+3.6*math.sin(a),y-8*math.sin(a)-3.6*math.cos(a)),(x-8*math.cos(a)-3.6*math.sin(a),y-8*math.sin(a)+3.6*math.cos(a))]
        p=c.beginPath();p.moveTo(pts[0][0],H-pts[0][1])
        for xx,yy in pts[1:]:p.lineTo(xx,H-yy)
        p.close();c.setFillColor(HexColor(color));c.drawPath(p,fill=1,stroke=0)
        svg.append(f'<polygon points="{" ".join(f"{xx},{yy}" for xx,yy in pts)}" fill="{color}"/>')
def box(x,y,w,h,title,rows,kind='prep',size=17):
    fill,border=colors[kind];c.setFillColor(HexColor(fill));c.setStrokeColor(HexColor(border));c.setLineWidth(1.5)
    c.roundRect(x,H-y-h,w,h,0,fill=1,stroke=1)
    svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="0" fill="{fill}" stroke="{border}" stroke-width="1.5"/>')
    assert stringWidth(title,"Arial-Bold",18)<w-24,(title,w)
    text(x+w/2,y+29,title,18,True)
    assert 55 + (len(rows)-1)*23 <= h-12,(title,h)
    for i,row in enumerate(rows):
        assert stringWidth(row,'Arial',size)<w-20,(row,w)
        text(x+w/2,y+55+i*23,row,size)


# Unified methodology figure. Model groups denote alternative specifications.
box(30,25,300,96,'Dengue and climate',[
'Weekly cases, rainfall,','temperature and humidity'],'data',18)
box(350,25,300,96,'Socioeconomic data',[
'Municipal development index','IDHM (2010)'],'data',18)
box(670,25,300,96,'Geography and mobility',[
'Municipal adjacency and regions','Distance and transport connectivity'],'data',17)
for x in [180,500,820]:line([(x,121),(x,143)])
line([(180,143),(820,143)]);line([(500,143),(500,166)],arrow=True)
box(150,166,700,116,'Data integration and feature construction',[
'92 municipalities in Rio de Janeiro; weekly observations, 2017-2023',
'Align records and address missing climate values',
'Construct lags; log-transform case predictors; standardize covariates'],'prep',18)
line([(500,282),(500,316)],arrow=True)
box(150,316,700,96,'Temporal evaluation design',[
'Training: 2017-2021; held-out testing: 2022-2023',
'Rolling-origin validation and climate-lag sensitivity analysis'],'prep',18)
line([(500,412),(500,440)],arrow=True)
box(150,440,700,76,'Bayesian model fitting (all candidate models)',[
'Negative-binomial likelihood; R-INLA estimation'],'model',18)
line([(500,516),(500,546)],arrow=True)
box(150, 466 + 80, 700, 96, 'Non-spatial model development (M0-M5)', ['Municipality and temporal effects; climate and IDHM', 'Lagged local cases; interpolation sensitivity'], 'model', 18)
line([(500, 562 + 80), (500, 596 + 80)], arrow=True)
box(150, 596 + 80, 700, 76, 'Spatial benchmarks (S1-S2)', ['BYM2 spatial effects and lagged cases in adjacent municipalities'], 'model', 18)
line([(500, 672 + 80), (500, 693 + 80)])
line([(180, 693 + 80), (820, 693 + 80)])
for x in [180, 500, 820]:
    line([(x, 693 + 80), (x, 718 + 80)], arrow=True)
box(30, 718 + 80, 300, 123, 'Spatial and mobility sensitivity', ['S3-S5', 'Distance-based neighboring cases', 'Road/fluvial and air mobility'], 'model', 17)
box(350, 718 + 80, 300, 123, 'Rainfall heterogeneity', ['S6, S8, S10, S11', 'Regional, municipal, temporal,', 'and additive municipal + temporal'], 'model', 16.5)
box(670, 718 + 80, 300, 123, 'Temperature heterogeneity', ['S7, S9', 'Regional and municipal', 'temperature effects'], 'model', 17)
for x in [180, 500, 820]:
    line([(x, 841 + 80), (x, 860 + 80)])
line([(180, 860 + 80), (820, 860 + 80)])
line([(500, 860 + 80), (500, 878 + 80)])
line([(257.5, 878 + 80), (742.5, 878 + 80)])
for x in [257.5, 742.5]:
    line([(x, 878 + 80), (x, 900 + 80)], arrow=True)
box(30, 900 + 80, 455, 99, 'Predictive assessment: all models', ['Held-out MAE, RMSE, WAPE and R-squared', 'Observed and predicted incidence'], 'eval', 18)
box(515, 900 + 80, 455, 99, 'Model assessment: all models', ['DIC and WAIC', "Residual spatial autocorrelation: Moran's I"], 'eval', 18)
for x in [257.5, 742.5]:
    line([(x, 999 + 80), (x, 1018 + 80)])
line([(257.5, 1018 + 80), (742.5, 1018 + 80)])
line([(500, 1018 + 80), (500, 1038 + 80)], arrow=True)
box(150, 1038 + 80, 700, 76, 'Rainfall inference', ['Posterior relative risks and credible intervals across space and time'], 'eval', 18)
c.showPage()
c.save()
svg.append('</svg>')
(out / 'dengue_study_workflow.svg').write_text('\n'.join(svg))
