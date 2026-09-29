from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from xml.sax.saxutils import escape
out='output/pdf/Ethan_Massey_Resume_Bank_of_America_SWE.pdf'
styles={
 'body':ParagraphStyle('body',fontName='Times-Roman',fontSize=10.2,leading=12.3),
 'bullet':ParagraphStyle('bullet',fontName='Times-Roman',fontSize=10.2,leading=12.3,leftIndent=10,firstLineIndent=-8,spaceAfter=2),
 'head':ParagraphStyle('head',fontName='Times-Bold',fontSize=11,leading=13,spaceBefore=8),
 'name':ParagraphStyle('name',fontName='Times-Bold',fontSize=19,leading=22,alignment=TA_CENTER),
 'contact':ParagraphStyle('contact',fontName='Times-Roman',fontSize=10,leading=13,alignment=TA_CENTER),
}
story=[]
def p(s,style='body'):return Paragraph(s,styles[style])
def add(s):story.append(p(s))
def section(s):
 story.append(p(s,'head'));story.append(HRFlowable(width='100%',thickness=.6,color=colors.black,spaceAfter=5))
def row(l,r):
 t=Table([[p(l),p(r)]],colWidths=[402,126]);t.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),0),('RIGHTPADDING',(0,0),(-1,-1),0),('TOPPADDING',(0,0),(-1,-1),0),('BOTTOMPADDING',(0,0),(-1,-1),0)]))
 t._cellvalues[0][1].style=ParagraphStyle('right',parent=styles['body'],alignment=2)
 story.append(t)
def bullet(s):story.append(p('&#8226; '+s,'bullet'))
story.append(p('Ethan Massey','name'))
story.append(p('Cambridge, Massachusetts | <link href="mailto:enm2375@mit.edu">enm2375@mit.edu</link> | <link href="http://www.linkedin.com/in/ethan-massey-b93812378">LinkedIn</link> | <link href="https://github.com/EthanMassey2007">GitHub</link>','contact'))
section('EDUCATION')
row('<b>Massachusetts Institute of Technology</b>','Cambridge, MA')
row('Bachelor of Science in Mathematics and Computer Science','Expected May 2028')
add('GPA: 4.8/5.0')
add('<b>Relevant Coursework:</b> Design and Analysis of Algorithms, Introduction to Algorithms, Fundamentals of Programming, Mathematics for Computer Science, Probability and Random Variables, Real Analysis')
story.append(Spacer(1,3))
row('<b>Dallas College</b>','Dallas, TX')
row('Associate of Science','May 2025')
add('GPA: 4.0/4.0; 80 Credit Hours')
section('TECHNICAL SKILLS')
add('<b>Programming Languages:</b> Python, Java')
add('<b>Software Development:</b> Git, Flask, Pygame, Pygbag, Folium; web applications and browser deployment')
add('<b>AI &amp; Data:</b> AI literacy fundamentals, prompt engineering, machine learning, data visualization, Bayesian modeling')
add('<b>Libraries &amp; Tools:</b> NumPy, pandas, scikit-learn, XGBoost, Optuna, SHAP')
add('<b>Certification:</b> Information Technology Specialist: Java, May 2024')
section('EXPERIENCE &amp; PROJECTS')
row('<b>MIT Senseable City Lab | Undergraduate Researcher</b>','Aug 2025 - Sep 2026')
add('<link href="https://github.com/EthanMassey2007/Hierarchical-Bayesian-Model">Research Repository</link>')
bullet('Built an end-to-end <b>Python data pipeline</b> to ingest, process, and integrate 15 years of epidemiological, climate, demographic, and geospatial data across <b>92 municipalities</b> for spatiotemporal forecasting.')
bullet('Developed an interactive <b>Flask/Folium web application</b> to visualize municipality-level forecasts, historical outbreaks, and spatial disease patterns.')
bullet('Implemented temporal-lag and graph-based spatial feature pipelines using municipality adjacency networks; automated model training and hyperparameter optimization with <b>XGBoost and Optuna</b>.')
bullet('Implemented hierarchical Bayesian models with MCMC and municipality and time-varying parameters, achieving <b>0.878 out-of-sample R<super>2</super></b> while quantifying uncertainty and changing spatiotemporal relationships.')
bullet('First author of a paper accepted to <b>IEEE IGARSS 2026</b> on spatiotemporal dengue forecasting; journal manuscript on spatiotemporal dengue modeling and Bayesian inference in preparation.')
story.append(Spacer(1,4))
row('<b>Graph Theory and Combinatorics Research</b>','Sep 2024 - Dec 2024')
bullet('Co-developed <b>Alytiqcon</b>, an interactive Python/Pygame graph-analysis application supporting construction and analysis of directed, undirected, weighted, and unweighted graphs.')
bullet('Implemented <b>Fleury\'s and Hierholzer\'s algorithms</b>, backtracking, and dynamic programming for graph traversal, pathfinding, and combinatorial search; <b>deployed the application to the browser</b> with Pygbag.')
bullet('Co-authored research on graph algorithms, combinatorial optimization, and real-world routing applications; invited to present at <b>TUMC 2024</b>.')
story.append(Spacer(1,4))
row('<b>M3 International Mathematics Competition</b>','Feb 2025')
bullet('Earned <b>Honorable Mention (top 20 of 800 submissions)</b> in the M3 Challenge; co-authored a 30-page modeling paper on heatwave risk, power demand, and infrastructure vulnerability in Memphis.')
bullet('Built Python models for indoor heat transfer, 50-year electricity-demand forecasting, and neighborhood vulnerability, including a neural network trained over <b>5,000 gradient-descent iterations</b> and sensitivity analysis.')
section('LEADERSHIP &amp; ACTIVITIES')
row('<b>MIT Native American and Indigenous Association | President</b>','Aug 2026 - Present')
bullet('Lead chapter programming, conference participation, industry outreach, and professional-development initiatives for MIT students.')
row('<b>AISES, MIT Chapter | Executive</b>','Aug 2026 - Present')
SimpleDocTemplate(out,pagesize=(612,792),rightMargin=42,leftMargin=42,topMargin=30,bottomMargin=30,title='Ethan Massey - Software Engineer Resume',author='Ethan Massey').build(story)
import pdfplumber
with pdfplumber.open(out) as pdf:
 print('Pages:',len(pdf.pages))
 for i,page in enumerate(pdf.pages):
  page.to_image(resolution=130).save(f'tmp/pdfs/resume-final-{i+1}.png')
 print(pdf.pages[0].extract_text())
