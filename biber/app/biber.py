import streamlit as st
import numpy as np
import pandas as pd
import pickle 
import io
import base64
import df_helper as helper # custom script 
import matplotlib.pyplot as plt
import seaborn as sns
from streamlit_tags import st_tags, st_tags_sidebar
import random
st.set_option('deprecation.showPyplotGlobalUse', False)
#matplotlib.use('TkAgg')
import sys
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, cut_tree
from scipy.spatial.distance import squareform
from matplotlib.pyplot import gcf
from fpdf import FPDF
import tempfile
import time
sys.setrecursionlimit(100000)
sns.set_palette("Paired") 
data = None


def highlight(txt):
    return '<span style="color: #F04E4E">{}</span>'.format(txt)


nat_sort = lambda l : ns.natsorted(l)
removeinvalidchars = lambda s : re.sub('[^a-zA-Z0-9\n\._ ]', '', s).strip()

#st.subheader("LOREM:")
#st.write('''* LOREM2.*''')


def download_file(df,name, extension):
    if extension == 'csv': # csv 
        csv = df.to_csv(index=True)
        b64 = base64.b64encode(csv.encode()).decode() 
    else: # pickle
        b = io.BytesIO()
        pickle.dump(df, b)
        b64 = base64.b64encode(b.getvalue()).decode()
    
    st.download_button('Download CSV of '+name, csv)

def create_download_link(val, filename):
    b64 = base64.b64encode(val)  # val looks like b'...'
    return f'<a href="data:application/octet-stream;base64,{b64.decode()}" download="{filename}.pdf">Download Plot</a>'


@st.cache()
def Genes(df,frac):
    if frac < 100:
        df = df.sample(frac=frac/100)
    return df
 

def genes_select(df,container):
   _col1, _col2, _colfoo = container.columns([1,1,3])
    
   frac=0.1
   
   _col1.header('Selected Genes')
   with _col1:
         statusSample = st.radio("", ('Paste a List','Random Sample'),1)
   
   with _col2:
    form = st.form(key='select_genes')
    if (statusSample  == 'Paste a List'):
         
         user_input = form.text_input(label='Space Separated Gene List',value="Bcin08g07060 Bcin01g07170 Bcin06g04710")
         user_input = user_input.split(" ")
         user_input = list(filter(None, user_input))
         if len(user_input) >= 2 or len(user_input)<=250:
            
            try:
               df=df.loc[user_input]
               form.form_submit_button(label='Submit')
                
            except KeyError:
               df = Genes(df,frac)
               st.button("The Selected Genes are not included in the experiments")
               
               form.form_submit_button(label='Submit')
         else:
            df = Genes(df,frac)
            form.form_submit_button(label='Submit')
        
    else:    # SAMPLE SIZE
    
        frac = form.slider('Random sample of genes (%)', 0.1, 2.0, frac)
        df = Genes(df,frac)
        form.form_submit_button(label='Submit')
   
   _colfoo.success(str(len(df))+' Genes Selected')
   genes = _colfoo.expander("Gene List ", expanded=False)
   genes.write(" ".join(df.index.to_list()))
   
   return df


@st.cache(allow_output_mutation=True)
def labels_palette(strains):
    cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
    # Label 1
     
    strain_labels = strains["botrytis_strain"]
    strain_pal = sns.cubehelix_palette(strain_labels.unique().size, light=.9, dark=.1, reverse=True, start=1, rot=-2)
    strain_lut = dict(zip(map(str, strain_labels.unique()), strain_pal))
    strain_colors = pd.Series(strain_labels, index=strains.index).map(strain_lut)


    node_labels = strains["culture_media"]
    node_pal = sns.color_palette('husl',node_labels.unique().size)
    node_lut = dict(zip(map(str, node_labels.unique()), node_pal))

    node_colors = pd.Series(node_labels, index=strains.index).map(node_lut)
    strain_node_colors = pd.DataFrame(strain_colors).join(pd.DataFrame(node_colors))
    

    plant_labels = strains["plant_material"]
    plant_pal = sns.color_palette("tab20c",plant_labels.unique().size)
    plant_lut = dict(zip(map(str, plant_labels.unique()), plant_pal))

    plant_colors = pd.Series(plant_labels, index=strains.index).map(plant_lut)
    strain_node_colors = pd.DataFrame(strain_node_colors).join(pd.DataFrame(plant_colors))
    return [strain_labels,strain_node_colors,strain_lut,node_labels,node_lut,plant_labels,plant_lut]


def heatmap(data,metadata,_container):
   ''' Plot HeatMap'''

   container = st.container()
   container.header('Heatmap set up')
   col1, col2,col3 = st.columns([0.3,0.3,0.3])
   with col1:
      form = st.form(key='Data Options')
      form.subheader('Data')
      log2 = form.radio("", ('Use Log2+1','No'),0)
      if log2 != 'No':
        data= np.log2(data+1) 
      
      form.subheader('Cluster the experiments?')
      clustExp = form.radio("", ('Yes','No'),1)
      col_cluster = False
      if clustExp != 'No':
           col_cluster = True
      
      form.subheader('Cluster the Genes?')
      k = form.slider('Clusters:', 2, 10, 4)
      form.form_submit_button(label='Apply')
   
   
    
   with col2:
      formSize = st.form(key='Options') 
      formSize.subheader('Size and color Scale')
      statusColor = formSize.radio("Color scheme: ", ('Continous', 'Categorized'),0)
      #width = formSize.slider("plot width", 15.1, 25., 25.)
      width=25.0
      height = formSize.slider("plot height", 16.0, 30., 16.)
   
      vmin = formSize.slider("min color range",1.0,data.values.max()/3.0,float(data.values.min()))
      vmax = formSize.slider("max color range",data.values.max()/3.0,data.values.max(),float(data.values.max()))  
      formSize.form_submit_button(label='Format')
   

   stride = 2
   lsize = 12
      
   # clustering
   Z = linkage(data, 'average')
   clusters = fcluster(Z, k, criterion='maxclust')
   
   # Create a categorical palette to identify the clusters
   clusters_pal = sns.husl_palette(len(clusters), s=.75)
   clusters_lut = dict(zip(clusters, clusters_pal))
    
   # Convert the palette to vectors that will be drawn on the side of the matrix
   clusters_colors = pd.Series(clusters,name='cluster').map(clusters_lut).to_numpy()
   # 
    #### COLOR PALETTE
   # conditional statement to print
   strain_labels,strain_node_colors,strain_lut,node_labels,node_lut,plant_labels,plant_lut=labels_palette(metadata)
   if (statusColor == 'Continous'):
      my_colors='YlGnBu'
   else:
      my_colors=[(0.2,0.3,0.3),(0.4,0.5,0.4),(0.1,0.7,0),(0.1,0.9,0)]
 
   kws = dict(square=True,cbar_kws=dict( orientation='horizontal'))   
   
   g=sns.clustermap(data, yticklabels=True, xticklabels=True, row_cluster=True,col_cluster=col_cluster, row_linkage=Z,col_linkage=None,vmin=vmin, vmax=vmax, col_colors=strain_node_colors,row_colors=clusters_colors,colors_ratio=0.01, cmap=my_colors , linewidth=0.005, linecolor='w',dendrogram_ratio=(.02, .2),figsize=(width, height), **kws )

   if g.dendrogram_row is None:
       print("Apparently, genes were not clustered.")
       return -1

   g.ax_row_dendrogram.remove() 
   
   g.ax_col_dendrogram.set_visible(False)
   g.ax_col_dendrogram.set_xlim([0,0])
   # Draw the legend bar for strain     
  
   legxx = []
   for label in node_labels.unique():
       x = g.ax_col_dendrogram.bar(0, 0, color=node_lut[label], label=label, linewidth=0)
       legxx.append(x)
   l2 = plt.legend(legxx, node_labels.unique(),ncol=2, loc="center", fontsize=13,handleheight=1.3,title='CULTURE MEDIA', bbox_to_anchor=(1.5, 1.5))
   
   legxx2 = []

   for label in plant_labels.unique():
       y = g.ax_col_dendrogram.bar(0, 0, color=plant_lut[label], label=label, linewidth=0)
       legxx2.append(y)
   l3 = plt.legend(legxx2, plant_labels.unique(), loc="center", ncol=2, handleheight=1.3, fontsize=13,title='PLANT MATERIAL', bbox_to_anchor=(0.6, 0.9), bbox_transform=gcf().transFigure)
   plt.gca().add_artist(l3)
   
   legxx3 = []

   for label in strain_labels.unique():
       z = g.ax_col_dendrogram.bar(0, 0, color=strain_lut[label], label=label, linewidth=0)
       legxx3.append(z)
   l4 = plt.legend(legxx3, strain_labels.unique(), loc="center", ncol=2, handleheight=1.3, fontsize=13,title='STRAIN', bbox_to_anchor=(0.8, 0.9), bbox_transform=gcf().transFigure)
   plt.gca().add_artist(l2)
   
   

   new_labels = []
   for l in g.ax_heatmap.axes.get_xticklabels():
      text= l.get_text()
      text = text+' -- '+metadata.loc[text,'experiment_description']
      l.set_text(text)
      new_labels.append(l)
   g.ax_heatmap.axes.set_xticklabels(new_labels,rotation = 90)
  
   
   ax = g.ax_heatmap
   x0, _y0, _w, _h = g.cbar_pos
   g.ax_cbar.set_position([x0, _y0+0.055, 0.25, 0.03])
   ax.set_xticklabels(ax.get_xmajorticklabels(),fontsize = lsize)
   g.ax_cbar.set_title('DeSeq2 UNITS')
   if (statusColor != 'Continous'):
      colorbar = ax.collections[0].colorbar
      colorbar.set_ticks([1.5 ,4.3,7.1 ,10.0 ])
      colorbar.set_ticklabels(['lower','low','high','higher'])
      g.ax_cbar.set_title('RELATIVE EXPRESSION')
   
   
   container4 = st.container()
   container4.pyplot(g)


   with col3:
        formSize = st.form(key='Download')
        formSize.subheader('Download')
        modeOptions = ['Show', 'Download PDF', 'View and Download Data']
        modeShow = formSize.radio("", modeOptions, index=0)
        formSize.form_submit_button(label='Submit')
        if modeShow=='Download PDF':
                
                fn = 'biber.pdf'+time.strftime("%Y%m%d-%H%M%S")
                img = io.BytesIO()
                g.savefig(img, format='pdf')
                btn = st.download_button(
                label="Download PDF",
                data=img,
                file_name=fn,
                mime="image/png"
                )
        


        elif modeShow=='View or Download Data':
                # reordering the index
                
                # heatmap indexs
                heatIndex = g.dendrogram_row.reordered_ind
                
                if clustExp != 'No':
                    heatColIndex = g.dendrogram_col.reordered_ind
                    data = data[data.columns[[heatColIndex]]]
                # reordre data index from up to bottom on heatmap
                data_reorder=data.iloc[heatIndex]
                # reorder clusters according to heatmap order
                clusters = [clusters[i] for i in heatIndex]
                data_reorder = pd.concat([pd.Series(clusters, index=data_reorder.index, name='Cluster'), data_reorder], axis=1)
                
                # reorder according original dataSet
                data_reorder=data_reorder.loc[data.index]
                st.write(data_reorder)
                download_file(data_reorder,'Genes' ,'csv')
                download_file(metadata,'Experiment Description','csv')
                return data_reorder
   return data


######
#####

def update_strain():
             st.session_state.all = 'all'
             st.session_state.idCounter+1


def increment_showFactor():
    st.session_state.showFactor = 1

st.session_state.all = 'ini'
st.session_state.all= 'all'
st.session_state.idCounter = 0
if 'showFactor' not in st.session_state:
    st.session_state.showFactor = 0



# outside class .. just to cache
@st.cache()
def factors():
   
   #metadata = pd.read_csv("/app/data/2022/metadata_beb_02082022.csv",sep=',')
   metadata = pd.read_csv("/app/data/production/metadata.csv",sep=',')
   metadata.set_index('experiment_ID', inplace=True)
   return metadata


class factor():
       zero = 0
       def __init__(self):
               self.factorName=None
               self.container=None
               self.text=None
               self.title=None
               self.df = factors()
               self.subfactors= None
               #self.idCounter = factor.idCounter
               #self.all = False

       def _multiselect(self):
          st.session_state.idCounter += 1

          _df=self.df[[self.factorName]]
          
          uniqueList = list(_df[self.factorName].unique())
          # check state
          st.session_state

          try:
            selected = st.multiselect(self.text, options=uniqueList, default=uniqueList, key=st.session_state.idCounter)
          except:
            selected = st.multiselect(self.text, options=st.session_state[st.session_state.idCounter], default=uniqueList, key=st.session_state.idCounter)
            

          self.df = self.df.loc[self.df[self.factorName].isin(selected)]
          if len(selected)<1:
             st.write(f'You must select at least one {self.factorName}.')
          else:
             self.zero = 1


       def set(self,factorName,text,title,all='j'):
            self.factorName=factorName
            self.text=text
            self.title=title
            self.all = all
            self._multiselect()


       def metadata(self):
           if self.zero == 0:
               metadata = self.df.sample(10,axis=0)

           else:
               metadata = self.df

           return self.df

       def __str__(self) -> str:
           return str(len(self.df))
      

def group_by(metadata):
    groupby_dict=metadata['group_for_averaging'].to_dict()
    return groupby_dict
	

def data_clean(df):
   ''' add anything here as run after Gene selection'''
   #st.success(str(len(df))+' Genes Selected')
   #genes = st.expander("See Genes data", expanded=False)
   #genes.table(df)
   return df



def flatten_dict(d):
    def items():
        for key, value in d.items():
            if isinstance(value, dict):
                for subkey, subvalue in flatten_dict(value).items():
                    yield key + "." + subkey, subvalue
            else:
                yield key, value[0]
    return dict(items())

def main():
    '''
	df <- experiments dataframe from metadata. Change on select
	metadata <- experiment descriptions
	'''
    
    st.set_page_config(page_title="BEB",layout='wide')
    #st.write(st.session_state)
    st.title('BEB')
    st.header('Botrytis Expression Browser')
    st.sidebar.image("/app/images/ibio.png", width=120) 
    st.write('<style>div.row-widget.stRadio > div{flex-direction:row;}</style>', unsafe_allow_html=True)    
    st.markdown(''' To visualize gene expression patterns in *Botrytis cinerea*, we developed BEB, a web-based *B. cinerea* gene Expression Browser. 
    BEB allows effortless visualization of transcript levels without advanced computational skills.
    To facilitate analysis, diverse clustering and visualization options can be applied. Please, go to the "Read Me" page (sidebar) for instructions and other options.''')
   
    st.markdown("""
        BEB is an Open Platform of the Millennium Institute of Integrative Biology [iBio][ibio] open-source and open science initiative.
                
        [ibio]: <https://www.ibio.cl/tecnologiaslibres/>
        """, unsafe_allow_html=True ) 

    st.markdown("""
        BEB gene-level count matrix (counts) and the metadata files can be found [here][here].
                
        [here]: <https://github.com/bioquimico/biber/tree/main/biber/app/data/production>
        """, unsafe_allow_html=True ) 

    # From GECO!
    st.sidebar.subheader('Mode')
    modeOptions = ['Analyze', 'Read Me']
    mode = st.sidebar.radio("", modeOptions, index=0)
    tabMethods = [mainPlot,readMe]
    tabMethods[modeOptions.index(mode)]()

def mainPlot():
    #st.markdown("<h2> BIBER </h2>", unsafe_allow_html=True)
    col1, mid, col2 = st.columns([1,1,20])
    with col1:
     pass 
     #   st.image('./images/change.jpg', width=64)
    with col2:
     #st.write('who')
     #st.image('./images/is.jpg', width=64)
     # html_string ='<h1><a href="http://cbv.unab.cl" > Paulo Canessa <small>Centro de Biotecnología Vegetal</small><a></h1>'
     # st.markdown(html_string, unsafe_allow_html=True)            
     pass
    container = st.sidebar.container()
    containerMid = st.container()
    containerMid2 = st.container()
    containerMid3 = st.container()
    
    experiments = factor()
    experiments.container= container
    metadata = experiments.metadata()




    try:
        
        
        with st.sidebar:
            with  st.form(key='selectmultiselect'):
             
                experiments.set('botrytis_strain','Strain','Strain') 
                experiments.set('tissue','Tissue','Tissue')
                
                if st.session_state.showFactor == 1: 
                    with st.expander("Other Factors", expanded=False):
                        experiments.set('culture_media','Culture Media','culture media')
                        experiments.set('treatment','Treatment','treatment')
                        experiments.set('plant_material','Plant Material','plant_material')
                        experiments.set('plant_tissue','Plant Tissue','plant tissue')
                        experiments.set('time','Time','time')
  
                 
                subm = st.form_submit_button(label='Select', on_click=increment_showFactor)
                
                                           
        
    except:
        container.success("No data to show with the selected restrains, please reload webpage")
    
    
    containerMid.header('Experiment Description') 
    col1,col2,col3 = containerMid.columns([1,1,1]) 
    col2.success(experiments.__str__()+'  Experimental replicates')

     
    dfData = helper.get_df()
    dfData = dfData[metadata.index.tolist()]
    dfData = genes_select(dfData,containerMid2)
     
    statusGroup = col1.radio("", ('Average of replicates','Individual replicates'),0)
    if statusGroup=='Individual replicates':
        statusGroup='No'

    if statusGroup != 'No':
       groupby_dict=metadata['group_for_averaging'].to_dict()
       groupby_dict2=metadata.groupby('group_for_averaging').apply(lambda g: g.index.tolist()).to_dict()
       dfData = dfData.groupby(groupby_dict, axis = 1).mean()
       dfData=pd.DataFrame(dfData).rename(columns=flatten_dict(groupby_dict2))
       col2.success(str(len(dfData.columns))+'  Experimental groups')
    
    experimentsExp = containerMid.expander("Experiments Metadata", expanded=False)
    experimentsExp.table(experiments.df)
    df_ord = heatmap(dfData,metadata,containerMid2)
    #download_file(df_ord,'Genes' ,'csv')
    #download_file(metadata,'Experiment Description','csv')


def readMe() :
    '''
    BEBS readme 
    '''
    
    #header.title('How to usage')
    st.markdown("""
        BEB is a [Streamlit] open-source app framework developed to visualize, explore, and detect gene expression patterns from multiple RNA-seq derived from curated *B. cinerea* experiments. For any questions, including database updates (see below), contact Paulo Canessa at [paulo.canessa@unab.cl](mailto:paulo.canessa@unab.cl) or Daniel Aguayo [bioquimico@gmail.com](mailto:bioquimico@gmail.com).
        
        [Streamlit]: <https://www.streamlit.io/>
        """, unsafe_allow_html=True )   

    st.markdown("""
        BEB have data for 11710 genes (e.g., Bcin02g02780) and 218 experiments (RNA-seq libraries, [preprint] version). Accordingly, experiments subsets could be selected by reducing the experimental factors (available in the left section of the screen), while genes by providing a list of "genes of interest (GOI)". For exploration or demonstrative purposes, BEB can also generate a random list of GOI. 
        [preprint]: <https://www.biorxiv.org/content/10.1101/2022.08.29.504976v1/>
        """, unsafe_allow_html=True )   
    
    st.markdown("""
        ### Quick Guide to Getting Started
        
        *** Factor selection ***

        1.	Select the *B. cinerea* Strains to be used for the analysis.
        1.	Select the Tissue to be used for the analysis.
        1.	The user can also select other factors for the analysis, such as culture media, treatment, plant material, and plant tissue.
        1.	By default, BEB displays the average of experimental replicates (when available).
        1.	It is worth noting that all factors are initially selected. Once a factor list has been defined, the experiments that fulfill the requirements will be displayed at the bottom.

        *** Gene selection ***
        
        1.	By default, BEB displays data for 12 random Genes, representing 0.01% of those available in the dataset. This quantity can be increased using the "Random sample of genes" slidebar.
        1.	Switch to BEB's Selected Gene mode by clicking "Paste a List."
        1.	Paste a space-separated Gene list, then submit. BEB utilizes the *B. cinerea* gene nomenclature, i.e. Bcin07g06770 Bcin04g00030 Bcin01g03930 Bcin10g02930 Bcin11g04790 Bcin05g08370.
        1.	Both experimental metadata and gene list will be displayed in the heatmap.

        *** Visualization Instructions ***

        1.	By default, BEB applies clusterization to all experimental conditions. The genes displaying similar expression patterns are denoted with a color code to the left of the heatmap. By default, BEB displays 4 clusters of genes, a number that can be modified using the "Cluster the Genes?" slidebar.
        1.	The user can also adjust the heatmap plot height and the min/max color range employed in the heatmap. To allow easy visualization of the expression data, the user can select two color schemes: the "Continuos" expression level (between min and max values) or the "Categorized", which automatically classifies all GOI expression values into quartiles. 
        1.	The "Download Plot" option activates a link to download the displayed heatmap, experimental descriptions and selected experiments data.
    
    """, unsafe_allow_html=True )

    st.markdown("""
        ### Database update and further (more complex) procedures. 
        BEB users are encouraged to reach out to us once they have uploaded their publication preprints and RNA-Seq expression data (e.g., SRA), such that we can upload their data to BEB.
        Shortly, underneath the online graphical user interface of BEB, there are two DBs that the server requires to generate each heatmap depicting gene expression: i) the one describing the experiments (straightforward; minutes to update); and ii) the more complex DB containing the read counts of each gene. 
        The latter requires more time to update. 
        
         *** Update Instructions: ***
        If any lab generates its RNA-seq data, there are two options that any user can employ to analyze this data:
        
        **1)** they can request us to map, quantify and upload the data to BEB's server
        
        **2)** they can execute a copy of BEB on a local computer system. If this option is chosen, the user needs to determine gene counts as indicated in the [preprint], add this piece of info to BEB's databases, and finally, explore the new data with the one pre-calculated on BEB. 
        Go to the [GitHub] information explaining how this task can be performed. 
        [preprint]: <https://www.biorxiv.org/content/10.1101/2022.08.29.504976v1/>
        [GitHub]: <https://github.com/ibioChile/CanessaLab/tree/master/Tools/>
        """, unsafe_allow_html=True )   



def genData():
    pass

def plotData():
	pass

main()
