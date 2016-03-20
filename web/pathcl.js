Ext.require(['*']);

Ext.onReady(function() {
   

    var cw;
    var currentPathCl=-1; // currently selected pathway cluster

    Ext.tip.QuickTipManager.init();    

    //Ext.state.Manager.setProvider(Ext.create('Ext.state.CookieProvider'));

    var clusterToolbar = Ext.create('Ext.toolbar.Toolbar', {
	    items: [ {
		    // xtype: 'button', // default for Toolbars
		    text: 'Show'
		},'->',{
		    xtype    : 'textfield',
		    icon: 'preview.png',
		    cls: 'x-btn-text-icon',
		    fieldLabel: 'number of genes',
		    labelStyle: 'white-space: nowrap;',
		    name     : 'nGenes',
		    emptyText: 'number of genes to show'
		}]
	});

    var tutorialWindow = new Ext.Window({
	id:'tutorial-window', 
	title: 'Video Tutorial',
	layout:'fit',  
	activeItem: 0,  

	defaults: {border:false}, 
	bbar: Ext.create('Ext.toolbar.Toolbar', {
	    padding: 5,
	    items   : [{
		xtype: 'checkbox',
		boxLabel: 'do not automatically show this tutorial video on startup when I visit next time',
		checked: Ext.util.Cookies.get("hidetutorial")!=null,
		name: 'dontshowtutorial',
		listeners: {
		    change: function(field, value) {
			if(value) {
			    var now = new Date();
			    var expiry = new Date(now.getTime() + 365 * 24 * 60 * 60 * 1000);
			    Ext.util.Cookies.set("hidetutorial",true,expiry)
			} else {
			    console.log("clearing hidetutorial");
			}
		    }
		}
	    },'->',{
		xtype: 'button',
		text: 'Close',
		handler: function() {
		    tutorialWindow.hide();
		}
	    }
           ]
	}),
	items : [{
		    id: "video",
		    html: '<iframe width="720" height="400" src="//www.youtube.com/embed/N2Ritx5yqrc?rel=0" frameborder="0" allowfullscreen></iframe>'
		}]	
	});

    

    /* DETAILED CLUSTERING VIEW */
    var detailMode=1; // 1: gene 2: pathway
    var detailItemList={}; // pathway or gene list
    var detailNGenes=20; // max genes to show 

    var detailPanel = Ext.create('Ext.panel.Panel', {
	bodyPadding: 5,
	autoScroll: true,
	showGenes: function(ids) {
	    if(ids === undefined || ids.length==0) return;
	    detailPanel.setLoading(true);
	    Ext.Ajax.request({
		url: 'genecl.json',
		method: 'POST',          
		waitTitle: 'Connecting',
		waitMsg: 'Sending data...',                                     
		params: {
		    "genes" : encodeURIComponent(JSON.stringify(ids))
		},
		scope:this,
		failure: function(r,o){console.log('failure:'); console.log(r);},
		success: function(response) {
		    var data = Ext.JSON.decode(response.responseText)
		    detailPanel.ids=ids;
		    detailPanel.genecldata=data;
		    detailPanel.mode=1;
		    Ext.getCmp('expressionDetailsPane').setTitle("Expression Details: selected genes");
		    detailPanel.redraw(detailPanel.genecldata)
		    detailPanelGearMenu.getComponent(0).disable()
		    //detailPanelGearMenu.getComponent(1).disable()
		    detailPanelGearMenu.getComponent(2).disable()
		    detailPanelGearMenu.getComponent(4).disable()
		    detailPanel.setLoading(false);
		}
	    });
	},
	showPathways: function(ids) {
	    if(ids === undefined || ids.length==0) return;
	    detailPanel.setLoading(true);
	    var ngenes=detailPanelGearMenu.getComponent(0).getValue();
	    var twosided=detailPanelGearMenu.getComponent(4).checked;
	    Ext.Ajax.request({
		url: 'pathwaygenes.json',
		method: 'POST',          
		waitTitle: 'Connecting',
		waitMsg: 'Sending data...',                                     
		params: {
		    "ngenes" : ngenes,
		    "twosided" : twosided,
		    "genes" : encodeURIComponent(JSON.stringify(ids)),
		    "trim" : hc.trim
		},
		scope:this,
		failure: function(r,o){console.log('failure:'); console.log(r);},
		success: function(response) {
		    var data = Ext.JSON.decode(response.responseText)
		    detailPanel.ids=ids;
		    detailPanel.genecldata=data;
		    detailPanel.mode=2; 
		    Ext.getCmp('expressionDetailsPane').setTitle("Expression Details: top genes in specified pathways");
		    detailPanel.redraw(detailPanel.genecldata)
		    detailPanelGearMenu.getComponent(0).enable()
		    //detailPanelGearMenu.getComponent(1).enable()
		    detailPanelGearMenu.getComponent(2).enable()
		    detailPanelGearMenu.getComponent(4).enable()
		    detailPanel.setLoading(false);
		}
	    });
	    
	},
	searchSimilar: function(pattern) { // request genes most closely matching current data.colcol
	    if(!detailPanel.hasOwnProperty('genecldata')) return;
	    detailPanel.setLoading(true);
	    if(pattern === undefined) {
		// determine the colcol or single gene
		if(detailPanel.genecldata.hasOwnProperty('colcols')) {
		    pattern=detailPanel.genecldata.colcols.data;
		} else {
		    pattern=detailPanel.genecldata.matrix.data.slice(0,detailPanel.genecldata.matrix.dim[1]);
		}
	    }
	    var ngenes=detailPanelGearMenu.getComponent(0).getValue();
	    var twosided=detailPanelGearMenu.getComponent(4).checked;
	    Ext.Ajax.request({
		url: 'patterngenes.json',
		method: 'POST',          
		waitTitle: 'Connecting',
		waitMsg: 'Sending data...',                                     
		params: {
		    "ngenes" : ngenes,
		    "twosided" : twosided,
		    "pattern" : encodeURIComponent(JSON.stringify(pattern)),
		    "trim" : hc.trim
		},
		scope:this,
		failure: function(r,o){console.log('failure:'); console.log(r);},
		success: function(response) {
		    var data = Ext.JSON.decode(response.responseText)
		    detailPanel.pattern=pattern;
		    detailPanel.genecldata=data;
		    detailPanel.mode=3;
		    Ext.getCmp('expressionDetailsPane').setTitle("Expression Details: genes matching specified pattern");
		    detailPanel.redraw(detailPanel.genecldata)
		    detailPanelGearMenu.getComponent(0).enable()
		    //detailPanelGearMenu.getComponent(1).enable()
		    detailPanelGearMenu.getComponent(2).enable()
		    detailPanelGearMenu.getComponent(4).enable()
		    detailPanel.setLoading(false);
		}
	    });
	},
	reload: function() { // goes back to the server to redraw the same ids
	    if(!detailPanel.hasOwnProperty('mode')) return;
	    switch(detailPanel.mode) {
	    case 1: detailPanel.showGenes(detailPanel.ids); break;
	    case 2: detailPanel.showPathways(detailPanel.ids);break;
	    case 3: detailPanel.searchSimilar(detailPanel.pattern);break;
	    default: console.log("reload requested with an undefined detailPanel mode");
	    }
	},
	redraw: function(data) { // redraws the panels without going back to the server
	    if(data === undefined) {
		if(detailPanel.hasOwnProperty('genecldata')) { data=detailPanel.genecldata; } else {  return; }
	    }
	    //$('.datapt').remove();
	    $('#geneclsvg .datapt').remove();
	    $('#geneclsvg').remove();
	    detailPanel.update("");    
	    s=detailPanel.getSize();
	    s.height=hc.geneUnitHeight*data.matrix.dim[0]+20; 
	    if(data.hasOwnProperty('colcols')) { s.height=s.height+hc.ccheight(data); }
	    var el = d3.select(detailPanel.getLayout().getElementTarget().dom)
	    var svg = el.append("svg").attr("id","geneclsvg").attr("width",(s.width-hc.padding.width)+"px").attr("height",(s.height-hc.padding.height)+"px").attr('xmlns','http://www.w3.org/2000/svg');
	    
	    if(data.hasOwnProperty('colcols')) {
		var colcol = colcolmap(svg.append("g").attr("transform","translate("+hc.hmleft()+","+hc.cctop(data)+")"), data.colcols, hc.hmwidth(s.width), hc.ccheight(data)); 
		colcol.append("title").text("1st principal component (PC) of the selected gene set expression: green - negative, white - neutral, orange - positive")
	    }
	    var cmap = colormap(svg.append("g").attr("id","gmapg").attr("transform","translate("+hc.hmleft()+","+hc.hmtop(data)+")"),data.matrix,hc.hmwidth(s.width),hc.hmheight(data,s.height)); 
	    if(data.hasOwnProperty('rowcols')) {
		var rowcols = sidecolormap(svg.append("g").attr("transform","translate("+hc.margins.left+","+hc.hmtop(data)+")"), data.rowcols, hc.rowcolWidth, hc.hmheight(data,s.height));
		rowcols.append("title").text("contribution of a gene to the 1st principal component: green - negative, white - neutral, orange - positive")
	    }
	    var rowLabelStep = hc.hmheight(data,s.height)/data.matrix.dim[0];
	    var rowlabg = svg.append("g")
		.attr("transform","translate("+hc.rowlableft(s.width)+","+hc.hmtop(data)+")")
		.append("g")
	    if(data.matrix.dim[0]==1) data.matrix.rows=[data.matrix.rows]

	    var rowlab = rowlabg
		.selectAll(".rowLabelg")
		.data(data.matrix.rows)
		.enter()
		.append("text")
		.text(function(d) { return(d); })
		.attr("x",0)
		.attr("y",function(d,i) { return Math.floor(i*rowLabelStep + rowLabelStep/2); })
		.classed("geneRowLabel",true)[0]
	    

	    var focus = svg.append("g")
		.attr("class","crosshair");
	    var focushl = focus.append("line").classed("crosshair",true).attr({"x1":hc.margins.left,"y1":Math.round(hc.hmtop(data)+cmap.y(1)/2),"x2":(hc.hmleft()+hc.hmwidth(s.width)),"y2":Math.round(hc.hmtop(data)+cmap.y(1)/2)});
	    var focusvl = focus.append("line").classed("crosshair",true).attr("id","genefocusvl").attr({"x1":Math.round(hc.hmleft()+cmap.x(1)/2),"y1":hc.cctop(data),"x2":Math.round(hc.hmleft()+cmap.x(1)/2),"y2":(hc.hmtop(data)+hc.hmheight(data,s.height))});
	    var focustx = focus.append("text").text("").attr("x",Math.round(hc.hmleft()+cmap.x(1)/2)).attr("y",Math.round(hc.hmtop(data)+cmap.y(1)/2));


	    cmap.g.append("rect").attr('x',0).attr('y',0).attr('width',hc.hmwidth(s.width)).attr('height',hc.hmheight(data,s.height)).attr('style','fill:white;fill-opacity:0.0;stroke:black;stroke-width:0.5;pointer-events:all;').attr("id","gmapev");
	    
	    var evg = d3.select("#gmapev");
	    evg.on("mousemove", function() {
		var coord=d3.mouse(evg.node());
		// figure out row and column index
		var bbox=evg.node().getBBox();
		var colIndex=Math.floor(coord[0]/bbox.width*data.matrix.dim[1])
		var rowIndex=Math.floor(coord[1]/bbox.height*data.matrix.dim[0])
		d3.select(rowlab[rowIndex]).classed('active', true);
		var newy=cmap.y(rowIndex); var newx=cmap.x(colIndex); 
		focushl.attr("transform","translate(0,"+newy+")");
		focusvl.attr("transform","translate("+newx+",0)");
		focustx.text("cell: "+data.matrix.cols[colIndex]);
		//focustx.attr("transform","translate("+coord[0]+","+(heatmapPos[1]+heatmapHeight-10)+")");
		if(newy>bbox.height/2) {  newy=newy-5; } else {   newy=newy+15; }
		if(newx>bbox.width/2) {
		    focustx.attr("transform","translate("+(newx-5)+","+newy+")");
		    focustx.attr("text-anchor","end")
		} else {
		    focustx.attr("transform","translate("+(newx+5)+","+newy+")");
		    focustx.attr("text-anchor","start")
		}
		d3.select('#pathclfocusvl').attr("transform","translate("+newx+",0)");
	    }).on("mouseenter",function() {
		el.classed('highlighting', true);
		d3.select('#pathclfocusvl').classed("highlighting",true)
	    }).on("mouseleave",function() {
		el.classed('highlighting', false);
		d3.select('#pathclfocusvl').classed("highlighting",false)
	    })

	},
	listeners: {
	    resize: function(cmp,width,height,oldWidth,oldHeight,opts) {
		cmp.redraw();
	    }
	}

    })


    // clear filter filed trigger button
    Ext.define('Ext.ux.CustomTrigger', {
	extend: 'Ext.form.field.Trigger',
	alias: 'widget.customtrigger',
	initComponent: function () {
            var me = this;
            me.triggerCls = 'x-form-clear-trigger';
            me.callParent(arguments);
	},
	// override onTriggerClick
	onTriggerClick: function() {
	    if(this.getValue()!='') {
		this.setRawValue('');
		this.fireEvent('change',this,'');
	    }
	}
    });


    /* PATHWAY CLUSTER INFO */

    Ext.define('clinfo',{
        extend: 'Ext.data.Model',
        fields: [
	    'name', 'id', {name: 'od',type: 'float'},{name: 'npc', type: 'integer'},{name: 'sign', type: 'integer'},{name: 'initsel', type: 'integer'}
        ],
	idProperty: 'id'
    });

    var clinfostore = Ext.create('Ext.data.Store', {
        id: 'clinfostore',
        model: 'clinfo',
        remoteSort: true,
        proxy: {
            type: 'jsonp',
            url: 'clinfo.json',
	    extraParams: {
		pathcl: currentPathCl
	    },
	    reader: {
		root: 'genes',
                totalProperty: 'totalCount'
            },
	    simpleSortMode: true,
        },
        autoLoad: false,
	pageSize: 50,
        remoteFilter: true,
	listeners: {
	    load: function(r) {
		// use the supplied 'initsel' to set the initial selection
		clSelectModel.suspendEvent('selectionchange')
		for(var i=0;i<r.data.items.length;i++) {
		    if(r.data.items[i].data.initsel==1) {
			clSelectModel.select(i,true);
		    }
		}
		clSelectModel.resumeEvent('selectionchange');
		clSelectModel.fireEvent('selectionchange',clSelectModel,clSelectModel.getSelection());
	    }
	}
    });

    var clSelectModel = Ext.create('Ext.selection.CheckboxModel', {
        listeners: {
            selectionchange: function(sm, selections) {
		var ids = $.map(selections,function(val,i) {
		    return(val.data.id)
		})
		if(ids.length>0) detailPanel.showPathways(ids);
            }
        }
    });
    
    var clInfoGrid = Ext.create('Ext.grid.Panel', {
        store: clinfostore,
	id: "clinfo",
	selModel: clSelectModel,
	height:'100%',
	columnLines:true,
	emptyText: 'No Matching Pathways',
        tbar: Ext.create('Ext.PagingToolbar', {
            store: clinfostore,
            displayInfo: false,
            //displayMsg: 'Displaying genes {0} - {1} of {2}',
            emptyMsg: "No pathways to display",
	    items:[
                {
		    flex:1,
		    width: 500,
		    minWidth: 50,
		    xtype: 'customtrigger',
		    emptyText: 'filter by name...',
		    listeners: {
			change: {buffer: 200, fn: function(field, value) {
			    if (value.length>0) {
				clinfostore.clearFilter(true);
				clinfostore.filter({property: 'name', value: value});
			    } else {
				clinfostore.clearFilter(false);
			    }
			}}
		    }
		}],
	    listeners: {
		afterrender: function() {
		    this.down('#refresh').hide();
		}
	    }
	}),
        columns: [
		  //{text: "name", flex: 1, dataIndex: 'name', sortable: true},
	    {
                text: 'overdispersion',
                dataIndex: 'od',
                flex:1,
                renderer: function (v, m, r) {
		    //m.tdAttr='data-qtip="'+r.data.name+'"';
                    var id = Ext.id();
                    Ext.defer(function () {
                        Ext.widget('progressbar', {
                            renderTo: id,
			    text: r.data.name,
                            value: v / 1,
                        });
                    }, 50);
		    if(r.data.sign=="1") {
			return Ext.String.format('<div class="positive" id="{0}"></div>', id);
		    } else {
			return Ext.String.format('<div class="negative" id="{0}"></div>', id);
		    }
                }
            },
            {text: "PC", width: 30, dataIndex: 'npc', sortable: true},
            /*{text: "overdispersion", width: 100, dataIndex: 'od', sortable: true}*/
        ]
    })



    /* GENE SET ENRICHMENT INFO */

    Ext.define('geinfo',{
        extend: 'Ext.data.Model',
        fields: [
	    'name', 'id', {name: 'fe',type: 'float'},{name: 'o', type: 'integer'},{name: 'u', type: 'integer'},{name: 'Z', type: 'float'},{name: 'Za', type: 'float'}
        ],
	idProperty: 'id'
    });

    var geinfostore = Ext.create('Ext.data.Store', {
        id: 'geinfostore',
        model: 'geinfo',
        remoteSort: true,
        proxy: {
            type: 'ajax',
            url: 'testenr.json',
	    actionMethods: {create: 'POST', read: 'POST', update: 'POST', destroy: 'POST'},
	    reader: {
		root: 'genes',
                totalProperty: 'totalCount'
            },
	    simpleSortMode: true,
        },
        autoLoad: false,
	pageSize: 50,
        remoteFilter: true,
    });

    var geSelectModel = Ext.create('Ext.selection.CheckboxModel', {
        listeners: {
            selectionchange: function(sm, selections) {
		var ids = $.map(selections,function(val,i) {
		    return(val.data.id)
		})
		if(ids.length>0) detailPanel.showPathways(ids);
            }
        }
    });
    
    var geInfoGrid = Ext.create('Ext.grid.Panel', {
        store: geinfostore,
	id: "geinfo",
	selModel: geSelectModel,
	height:'100%',
	columnLines:true,
	emptyText: 'No Enriched Pathways',
        tbar: Ext.create('Ext.PagingToolbar', {
            store: geinfostore,
            displayInfo: false,
            //displayMsg: 'Displaying genes {0} - {1} of {2}',
            emptyMsg: "No pathways to display",
	    items:[
                {
		    flex:1,
		    width: 500,
		    minWidth: 50,
		    xtype: 'customtrigger',
		    emptyText: 'filter by name...',
		    listeners: {
			change: {buffer: 200, fn: function(field, value) {
			    if (value.length>0) {
				geinfostore.clearFilter(true);
				geinfostore.filter({property: 'name', value: value});
			    } else {
				geinfostore.clearFilter(false);
			    }
			}}
		    }
		}],
	    listeners: {
		afterrender: function() {
		    this.down('#refresh').hide();
		}
	    }
	}),
        columns: [
            {text: "Pathway", flex: 1, dataIndex: 'name', sortable: true },
	    {text: "FE", width: 60, dataIndex: 'fe', sortable: true, tooltip: "fold enrichment"},
            {text: "Z", width: 60, dataIndex: 'Z', sortable: true, tooltip: "enrichment Z-score"},
            {text: "cZ", width: 60, dataIndex: 'Za', sortable: true, tooltip: "enrichment Z-score, corrected for multiple hypothesis testing"},
            {text: "n", width: 50, dataIndex: 'n', sortable: true, hidden: true, tooltip: "number of genes found in this pathway"},
	    {text: "u", width: 50, dataIndex: 'u', sortable: true, hidden: true, tooltip: "total number of genes annotated for this pathway"}
        ]
    })


    /* gene info tab */
    Ext.define('ginfo',{
        extend: 'Ext.data.Model',
        fields: [
	  'gene', {name: 'var',type: 'float'},{name: 'svar', type: 'float'}
        ],
	idProperty: 'gene'
    });

    var ginfostore = Ext.create('Ext.data.Store', {
        id: 'ginfostore',
        model: 'ginfo',
        remoteSort: true,
        proxy: {
            type: 'jsonp',
            url: 'genes.json',
	    reader: {
		root: 'genes',
                totalProperty: 'totalCount'
            },
	    simpleSortMode: true,
        },
        sorters: [{
            property: 'var',
            direction: 'DESC'
        }],
	pageSize: 100,
        remoteFilter: true,
        autoLoad: true
    });



    var geneSelModel = Ext.create('Ext.selection.CheckboxModel', {
        listeners: {
            selectionchange: function(sm, selections) {
		    var ids = $.map(selections,function(val,i) {
			    return(val.data.gene)
			})
		    detailPanel.showGenes(ids);
		    //grid4.down('#removeButton').setDisabled(selections.length === 0);
            }
        }
    });

    var geneGrid = Ext.create('Ext.grid.Panel', {
        store: ginfostore,
	selModel: geneSelModel,
        columns: [
            {text: "Gene", flex: 1, dataIndex: 'gene', sortable: true,
	     renderer: function(value) {
		 return Ext.String.format('<a href="http://www.informatics.jax.org/searchtool/Search.do?query={0}" target="_blank">{1}</a>',value,value)
	     }
	    },
            {text: "Variance", width: 80, dataIndex: 'var', sortable: true}
        ],
	//features: [filters],
	height:'100%',
	columnLines:true,
	emptyText: 'No Matching Genes',
	//forceFit: true
        //renderTo:'example-grid',
        //width: 800,
        //height: 300
	// paging bar on the bottom
        tbar: Ext.create('Ext.PagingToolbar', {
            store: ginfostore,
            displayInfo: false,
            //displayMsg: 'Displaying genes {0} - {1} of {2}',
            emptyMsg: "No genes to display",
	    items:[
                {
		    flex:1,
		    width: 500,
		    minWidth: 50,
		    xtype: 'customtrigger',
		    emptyText: 'filter by gene name...',
		    listeners: {
			change: {buffer: 200, fn: function(field, value) {
			    if (value.length>0) {
				ginfostore.clearFilter(true);
				ginfostore.filter({property: 'gene', value: value});
			    } else {
				ginfostore.clearFilter(false);
			    }
			}}
		    }
		}],
	    listeners: {
		afterrender: function() {
		    this.down('#refresh').hide();
		}
	    }
	}),
	listeners: {
	    viewready: function() {
		    // select top 20 genes 
		    this.selModel.suspendEvent('selectionchange')
		    for(var i=0;i<Math.min(this.store.data.items.length,20);i++) {
			this.selModel.select(i,true);
		    }
		    this.selModel.resumeEvent('selectionchange');
		    this.selModel.fireEvent('selectionchange',this,this.selModel.getSelection());
	    }
	}
    });



    /* gene info tab */
    Ext.define('pinfo',{
        extend: 'Ext.data.Model',
        fields: [
	  'id','name', {name: 'Z',type: 'float'},{name: 'aZ',type: 'float'},{name: 'score',type: 'float'},{name: 'n', type: 'integer'},{name: 'npc', type: 'integer'}
        ],
	idProperty: 'id'
    });

    var pinfostore = Ext.create('Ext.data.Store', {
        id: 'pinfostore',
        model: 'pinfo',
        remoteSort: true,
        proxy: {
            type: 'jsonp',
            url: 'pathways.json',
	    reader: {
		root: 'genes',
                totalProperty: 'totalCount'
            },
	    simpleSortMode: true,
        },
	pageSize: 100,
        remoteFilter: true,
        autoLoad: true
    });

    var pathwaySelModel = Ext.create('Ext.selection.CheckboxModel', {
        listeners: {
            selectionchange: function(sm, selections) {
		var ids = $.map(selections,function(val,i) {
		    return(val.data.id)
		})
		detailPanel.showPathways(ids);
            }
        }
    });

    var pathwayGrid = Ext.create('Ext.grid.Panel', {
        store: pinfostore,
	selModel: pathwaySelModel,
        columns: [
            {text: "Pathway", flex: 1, dataIndex: 'name', sortable: true, tooltip: "pathway / gene set name" },
            {text: "Z", width: 60, dataIndex: 'Z', sortable: true, hidden: true, tooltip: "overdispersion Z score"},
            {text: "cZ", width: 60, dataIndex: 'aZ', sortable: true, tooltip: "overdispersion Z score, adjusted for multiple hypothesis"},
	    {text: "co.Z", width: 60, dataIndex: 'sh.Z', sortable: true, hidden: true, tooltip: "pathway coherence Z-score"},
            {text: "co.cZ", width: 60, dataIndex: 'sh.aZ', sortable: true, hidden: true, tooltip: "pathway coherence Z-score, corrected for multiple hypothesis testing"},
            {text: "n", width: 60, dataIndex: 'n', sortable: true, hidden: true, tooltip: "number of genes in the pathway"},
	    {text: "nPC", width: 60, dataIndex: 'npc', sortable: true, hidden: true, tooltip: "principal component number"},
	    {text: "score", width: 60, dataIndex: 'score', sortable: true, tooltip: "observed/expected overdispersion"}
        ],
	//features: [filters],
	emptyText: 'No Matching Pathways',
        height: '100%',
	columnLines:true,
	// paging bar on the bottom
        tbar: Ext.create('Ext.PagingToolbar', {
            store: pinfostore,
            displayInfo: false,
            //displayMsg: 'Displaying genes {0} - {1} of {2}',
            emptyMsg: "No pathways to display",
	    items:[
                {
		    flex:1,
		    width: 500,
		    minWidth: 50,
		    xtype: 'customtrigger',
		    emptyText: 'filter by pathway name...',
		    listeners: {
			change: {buffer: 200, fn: function(field, value) {
			    if (value.length>0) {
				pinfostore.clearFilter(true);
				pinfostore.filter({property: 'name', value: value});
			    } else {
				pinfostore.clearFilter(false);
			    }
			}}
		    }
		}],
	    listeners: {
		afterrender: function() {
		    this.down('#refresh').hide();
		}
	    }
	})
    });



    var infotab = Ext.create('Ext.tab.Panel', {
	    //tabPosition: 'right',
	    defaults: {
		bodyPadding: 0,
		layout: 'fit',
		iconCls: 'tab-icon'
	    },
	    items: [
		{ title: 'Pathways', items:[pathwayGrid] },
		{ title: 'Genes', items:[geneGrid] },
		{ title: 'Cluster', items:[clInfoGrid], itemId: 'clustertab', hidden:true },
		{ title: 'Enrichment', items:[geInfoGrid], itemId: 'enrichmenttab', hidden:true }
	    ]
	});

    
    function colcolmap(svg, data, width, height) {
	// Check for no data
	if (data.length === 0)
	    return function() {};
	var cols = data.dim[1]
	var rows = data.dim[0];
	
	var merged = data.data;
	
	var x = d3.scale.linear().domain([0, cols]).range([0, width]);
	var y = d3.scale.linear().domain([0, rows]).range([0, height]);

	var g = svg.append("g").classed("colormap",true);
	var rect = g.selectAll("rect").data(merged);
	rect.enter().append("rect").classed("colcoldatapt", true);
	//rect.exit().remove();
	var cccmap= function(d) { return d; };
	if(data.hasOwnProperty('domain')) {
	    var color = d3.scale.linear()
		.domain(data.domain)
		.range(data.colors);
	    cccmap= function(d) { return color(d); };
	};

	rect.attr("x", function(d, i) {
		return x(i % cols);
	    })
	    .attr("y", function(d, i) {
		return y(Math.floor(i / cols));
	    })
	    .attr("width", x(1))
	    .attr("height", y(1))
	    .attr("fill", cccmap);
	svg.append("rect").attr('x',0).attr('y',0).attr('width',width).attr('height',height).attr('style','fill:none;stroke:black;stroke-width:0.5;');
	return g;
    };
    function sidecolormap(svg, data, width, height) {
	// Check for no data
	if(data === undefined) 
	    return function() {};
	    
	var cols = data.dim[1]
	var rows = data.dim[0];
	
	var merged = data.data;
	
	var x = d3.scale.linear().domain([0, cols]).range([0, width]);
	var y = d3.scale.linear().domain([0, rows]).range([0, height]);
	var color = d3.scale.linear()
	    .domain(data.domain)
	    .range(data.colors);
	

	var g = svg.append("g").classed("colormap",true);
	var rect = g.selectAll("rect").data(merged);
	rect.enter().append("rect").classed("datapt", true);
	//rect.exit().remove();
	rect.property("colIndex", function(d, i) { return i % cols; })
	    .property("rowIndex", function(d, i) { return Math.floor(i / cols); })
	    .attr("x", function(d, i) {
		return x(i % cols);
	    })
	    .attr("y", function(d, i) {
		return y(Math.floor(i / cols));
	    })
	    .attr("width", x(1))
	    .attr("height", y(1))
	    .attr("fill", function(d) { return color(d); });
	    //.append("title").text(function(d) { return d + ""; });
	svg.append("rect").attr('x',0).attr('y',0).attr('width',width).attr('height',height).attr('style','fill:none;stroke:black;stroke-width:0.5;');
	return g;
    };

    function colormap(svg, data, width, height) {
	// Check for no data
	if (data.length === 0)
	    return function() {};
	var cols = data.dim[1]
	var rows = data.dim[0];
	var merged = data.data;
	
	var x = d3.scale.linear().domain([0, cols]).range([0, width]);
	var y = d3.scale.linear().domain([0, rows]).range([0, height]);
	var color = d3.scale.linear()
	    .domain(data.domain)
	    .range(data.colors);
	
	
	var g = svg.append("g").classed("colormap",true);
	//g.append("rect").attr('x',0).attr('y',0).attr('width',width).attr('height',height).attr('style','fill:white;stroke:none;');
	//svg.append("rect").attr('x',0).attr('y',0).attr('width',width).attr('height',height).attr('style','fill-opacity:1;fill:black;stroke:none;pointer-events:all;').classed("mapcatcher",true);
	var rect = g.selectAll("rect").data(merged);
	rect.enter().append("rect").classed("datapt", true);
	//rect.exit().remove();
	rect.attr("x", function(d, i) {
		return x(i % cols);
	    })
	    .attr("y", function(d, i) {
		return y(Math.floor(i / cols));
	    })
	    .attr("width", x(1))
	    .attr("height", y(1))
	    .attr("fill", function(d) { return color(d); });
	    //.append("title").text(function(d) { return d + ""; });
	svg.append("rect").attr('x',0).attr('y',0).attr('width',width).attr('height',height).attr('style','fill:none;stroke:black;stroke-width:0.5;');
	return {
	    x: x,
	    y: y,
	    g: g
	};
    };



    function updatePathclInfo(pathcl) {
	if(pathcl == currentPathCl) return;
	clinfostore.getProxy().setExtraParam("pathcl",pathcl)
	clinfostore.load();
	infotab.child("#clustertab").tab.show();
	infotab.setActiveTab(2);
	infotab.getActiveTab().setTitle("Cluster "+pathcl)
	currentPathCl=pathcl;
    }
    
    /* heatmap config */
    var hc = {
	spacing: 2,
	geneUnitHeight: 15,
	margins: {top:2,right:60,bottom:2,left:1},
	padding: {width:28,height:10},
	colcolUnitHeight: 10,
	rowcolWidth: 10,
	colDendHeight: 50,
	trim: 0,
	// heatmap left position
	hmleft: function() { return(this.margins.left+this.rowcolWidth+this.spacing)},
	// colcol top position
	cctop: function(data) { return(data.hasOwnProperty('coldend')? this.margins.top+this.colDendHeight+this.spacing : this.margins.top) },
	// colcol height
	ccheight: function(data) { return(this.colcolUnitHeight*data.colcols.dim[0]) },
	// heatmap top position
	hmtop: function(data) { 
	    if(data.hasOwnProperty('colcols')) {
		return(data.hasOwnProperty('coldend')? this.margins.top+this.colDendHeight+this.spacing*2 + this.ccheight(data) : this.margins.top+this.spacing*2 + this.ccheight(data))
	    } else {
		return(data.hasOwnProperty('coldend')? this.margins.top+this.colDendHeight+this.spacing : this.margins.top)
	    }
	},
	// heatmap hight
	hmheight: function(data,totalHeight) {
	    return((totalHeight-this.padding.height*2) - this.hmtop(data));
	},
	// heatmap width
	hmwidth: function(totalWidth) { return(totalWidth-this.margins.left-this.rowcolWidth-this.spacing-this.margins.right-this.padding.width*2)},
	rowlableft: function(totalWidth) {
	    return(this.hmleft()+this.hmwidth(totalWidth)+this.spacing);
	},
	rowlabelsize: function(data,totalHeight) {
	    return(Math.round(this.hmheight(data,totalHeight)/data.matrix.dim[0]*0.97*72/96))
	    //return(18);
	}
    };

    var clusterPanel = Ext.create('Ext.panel.Panel', {
	layout: 'fit',
	bodyPadding: 5,
	//ohtml: "Pathway Clustering",
	reload: function() {
	    clusterPanel.setLoading(true);
	    d3.json('pathcl.json',function(error,data) {
		if(error) { console.log(error); return;}
		clusterPanel.pathcldata=data;
		clusterPanelGearMenu.getComponent(0).suspendEvents();
		clusterPanelGearMenu.getComponent(0).setValue(data.matrix.domain[data.matrix.domain.length-1]);
		clusterPanelGearMenu.getComponent(0).setMaxValue(data.matrix.range[1]);
		clusterPanelGearMenu.getComponent(0).resumeEvents()
		if(data.hasOwnProperty('trim')) {
		    hc.trim=data.trim;
		    detailPanelGearMenu.getComponent(2).suspendEvents();
		    detailPanelGearMenu.getComponent(2).setValue(data.trim);
		    detailPanelGearMenu.getComponent(2).resumeEvents();
		}
		clusterPanel.redraw(clusterPanel.pathcldata);
		clusterPanel.setLoading(false);
	    })
	},
	redraw: function(data) {
	    if(data === undefined) {
		if(clusterPanel.hasOwnProperty('pathcldata')) { data=clusterPanel.pathcldata; } else { clusterPanel.reload(); return; }
	    }
	    $('#pathclsvg .datapt').remove();
	    $('#pathclsvg').remove();
	    clusterPanel.update("");
	    s=clusterPanel.getSize();
	    var el = d3.select(clusterPanel.getLayout().getElementTarget().dom)
	    var svg = el.append("svg").attr("id","pathclsvg").attr("width",(s.width-hc.padding.width)+"px").attr("height",(s.height-hc.padding.height)+"px").attr('xmlns','http://www.w3.org/2000/svg');
	    
	    var dg = svg.append("g").attr("id","coldend").attr("transform", "translate(" + hc.hmleft() + " " + hc.margins.top + ") " +"scale(" + (hc.hmwidth(s.width)/72) + " " + (hc.colDendHeight/72) + ")");
	    // a workaround to append SVG elements into DOM
	    $(el.node()).append('<svg id="dummy" style="display:none"><defs>' + data.coldend + '</defs></svg>');
	    $(dg.node()).append($("#dummy g"));
	    $("#dummy").remove();
	    // colcol
	    var colcol = colcolmap(svg.append("g").attr("transform","translate("+hc.hmleft()+","+hc.cctop(data)+")"), data.colcols, hc.hmwidth(s.width), hc.ccheight(data)); 
	    colcol.append("title").text("custom cell classification colors")
	    // main heatmap
	    var cmap = colormap(svg.append("g").attr("id","cmapg").attr("transform","translate("+hc.hmleft()+","+hc.hmtop(data)+")"),data.matrix,hc.hmwidth(s.width),hc.hmheight(data,s.height)); 
	    // side colors
	    var rowcols = sidecolormap(svg.append("g").attr("transform","translate("+hc.margins.left+","+hc.hmtop(data)+")"), data.rowcols, hc.rowcolWidth, hc.hmheight(data,s.height));
	    rowcols.append("title").text("Overdispersion: white - low, black - high")
	    
	    // row labels
	    var rowlabelsize=hc.rowlabelsize(data,s.height);
	    var rowLabelStep = hc.hmheight(data,s.height)/data.matrix.dim[0];
	    var rowlabg = svg.append("g")
		.attr("transform","translate("+hc.rowlableft(s.width)+","+hc.hmtop(data)+")")
		.append("g")
	    var rowlab = rowlabg
		.selectAll(".rowLabelg")
		.data(data.matrix.rows)
		.enter()
		.append("text")
		.text(function(d) { return(d); })
		.style("font-size",rowlabelsize+"pt")
		.attr("x",0)
		.attr("y",function(d,i) { return Math.floor(i*rowLabelStep + rowLabelStep/2); })
		.classed("rowLabel",true)[0]
	    
	    var focus = svg.append("g").attr("class","crosshair");
	    var focushl = focus.append("line").classed("crosshair",true).attr({"x1":hc.margins.left,"y1":Math.round(hc.hmtop(data)+cmap.y(1)/2),"x2":(hc.hmleft()+hc.hmwidth(s.width)),"y2":Math.round(hc.hmtop(data)+cmap.y(1)/2)});
	    var focusvl = focus.append("line").classed("crosshair",true).attr("id","pathclfocusvl").attr({"x1":Math.round(hc.hmleft()+cmap.x(1)/2),"y1":hc.cctop(data),"x2":Math.round(hc.hmleft()+cmap.x(1)/2),"y2":(hc.hmtop(data)+hc.hmheight(data,s.height))});
	    var focustx = focus.attr("id","pathclfocustx").append("text").text("").attr("x",Math.round(hc.hmleft()+cmap.x(1)/2)).attr("y",Math.round(hc.hmtop(data)+cmap.y(1)/2));

	    cmap.g.append("rect").attr('x',0).attr('y',0).attr('width',hc.hmwidth(s.width)).attr('height',hc.hmheight(data,s.height)).attr('style','fill:white;fill-opacity:0.0;stroke:black;stroke-width:0.5;pointer-events:all;').attr("id","cmapev");
	    
	    var evg = d3.select("#cmapev");
	    evg.on("click", function() {
		var coord=d3.mouse(evg.node());
		var bbox=evg.node().getBBox();
		var colIndex=Math.floor(coord[0]/bbox.width*data.matrix.dim[1])
		var rowIndex=Math.floor(coord[1]/bbox.height*data.matrix.dim[0])
		updatePathclInfo(rowlab[rowIndex].textContent)
	    }).on("mousemove", function() {
		var coord=d3.mouse(evg.node());
		// figure out row and column index
		var bbox=evg.node().getBBox();
		var colIndex=Math.floor(coord[0]/bbox.width*data.matrix.dim[1])
		var rowIndex=Math.floor(coord[1]/bbox.height*data.matrix.dim[0])
		d3.select(rowlab[rowIndex]).classed('active', true);
		var newy=cmap.y(rowIndex); var newx=cmap.x(colIndex); 
		focushl.attr("transform","translate(0,"+newy+")");
		focusvl.attr("transform","translate("+newx+",0)");
		focustx.text("cell: "+data.matrix.cols[colIndex]);
		//focustx.attr("transform","translate("+coord[0]+","+(heatmapPos[1]+heatmapHeight-10)+")");
		if(newy>bbox.height/2) {  newy=newy-5; } else {   newy=newy+15; }
		if(newx>bbox.width/2) {
		    focustx.attr("transform","translate("+(newx-5)+","+newy+")");
		    focustx.attr("text-anchor","end")
		} else {
		    focustx.attr("transform","translate("+(newx+5)+","+newy+")");
		    focustx.attr("text-anchor","start")
		}
		d3.select('#genefocusvl').attr("transform","translate("+newx+",0)");
	    }).on("mouseenter",function() {
		el.classed('highlighting', true);
		d3.select('#genefocusvl').classed("highlighting",true)
	    }).on("mouseleave",function() {
		el.classed('highlighting', false);
		d3.select('#genefocusvl').classed("highlighting",false)
	    })
	    
	},
	listeners: {
	    resize: function(cmp,width,height,oldWidth,oldHeight,opts) {
		cmp.redraw(cmp.pathcldata);
	    }
	}
    });


    var clusterPanelGearMenu = Ext.create('Ext.menu.Menu', {
        id: 'clusterGearMenu',
        style: {
            overflow: 'visible'     // For the Combo popup
        },
        items: [{
		fieldLabel: 'Z limit',
		name: 'zlim',
		xtype: 'numberfield',
	        value: -1,
		decimalPrecision: 3,
		minValue: 0.0,
	        maxValue: 100,
		width: 200,
	        disabled: false,
		tooltip: 'Set the range of overdispersion scores illustrated by colors',
		listeners : {
		    change : {buffer: 800, fn:function(f,v) {
			clusterPanel.pathcldata.matrix.domain=$.map($(Array(clusterPanel.pathcldata.matrix.domain.length)),function(val, i) { return (2*i*v/clusterPanel.pathcldata.matrix.domain.length - v); })
			clusterPanel.redraw()
		    }}
		}
	}
	],
	listeners:{
	    'mouseleave': {buffer: 1000, fn:function( menu, e, eOpts){
		menu.hide();
	    }}
	}
    });

    
    var detailPanelGearMenu = Ext.create('Ext.menu.Menu', {
        id: 'detailGearMenu',
        style: {
            overflow: 'visible'     // For the Combo popup
        },
        items: [{
		fieldLabel: 'N genes',
		xtype: 'numberfield',
	        tooltip: 'Number of genes to show in the Expression Details panel',
		label: 'N genes',
		value: 20,
		minValue: 1,
		maxValue: 1000,
		disabled: true,
		listeners : {
		    change : {buffer: 800, fn:function(f,v) {detailPanel.reload()}}
		}
/*	    }, '-',{
		text: 'Color Z-limit:',
		tooltip: 'Maximum Z score to determine the color range',
		canActivate:false

	    },{
		xtype: 'slider',
		label: 'N genes',
		value: 0.5,
		increment: 1,
		minValue: 0,
		maxValue: 100,
		tipText: function(thumb){ return Ext.String.format('{0}', (thumb.value/100*3.6).toFixed(2)); },
*/
	    }, {
		fieldLabel: 'Row height',
		name: 'rowheight',
		xtype: 'numberfield',
		value: hc.geneUnitHeight,
		minValue: 3,
		width: 200,
		maxValue: 100,
		tooltip: 'Number of pixels to use for each gene row',
		listeners : {
		    change : {buffer: 800, fn:function(f,v) {hc.geneUnitHeight=v; detailPanel.redraw()}}
		}
	    },{
		fieldLabel: 'Trim',
		name: 'trim',
		xtype: 'numberfield',
		value: hc.trim,
		decimalPrecision: 5,
		minValue: 0.0,
		maxValue: 0.5,
		width: 200,
		maxValue: 100,
		disabled: true,
		tooltip: 'Winsorization trim fraction',
		listeners : {
		    change : {buffer: 800, fn:function(f,v) {hc.trim=v; detailPanel.reload()}}
		}
            }, '-',
	    {
                text: 'High/low genes',
                checked: false,
		tooltip: 'Whether to include genes from both sides of the PC loading (true) or just high magnitude  (false)',
		disabled: true,
		listeners : {
		    checkchange : function() {detailPanel.reload()}
		}
	    }],
	listeners:{
	    'mouseleave': {buffer: 1000, fn:function( menu, e, eOpts){
		menu.hide();
	    }}
 }
    });

    var ngenesSlider = Ext.create('Ext.slider.Single', {
	label: 'N genes',
	tip: 'number of genes to show',
	tipText: function(thumb){ return Ext.String.format('show {0} genes', thumb.value); },
	width: 100,
	value: 20,
	increment: 1,
	minValue: 0,
	maxValue: 500,

    });

    var viewport = Ext.create('Ext.Viewport', {
        layout: {
            type: 'border',
            padding: 5
        },
        defaults: {
            split: true
        },
        items: [{ region: 'center',
		  layout: 'border',
		  items: [{ region: 'north',
			    layout: 'fit',
			    id: 'clusterPanel',
			    title: 'Pathway Clustering',
			    tools: [
				{ type:'gear',
				  tooltip: 'Settings',
				  handler: function(e, el,o,t) {
				      clusterPanelGearMenu.showBy(t);
				  }
				},
				{ type:'help',
				  tooltip: 'Tutorial',
				  handler: function(e, el,o,t) {
				      tutorialWindow.show(); 
				  }
				},
				{ type:'save',
				  tooltip: 'Save image as SVG file',
				  handler: function(e,el,o,t) {
				      var svg=d3.select("#pathclsvg");
				      if(!svg.empty()) {
					  // update some visual attributes
					  svg.selectAll("#coldend * path").style("vector-effect","non-scaling-stroke");
					  svg.selectAll(".rowLabel").style("dominant-baseline","middle").style("font","18px sans-serif");
					  svg.select("#pathclfocustx").text("");
					  
					  var b64 = window.btoa(svg.node().parentNode.innerHTML);
					  writeAndClickLink('data:application/octet-stream;base64,\n'+b64,'pathway_clusters.svg')
				      }
				  }
				}
			    ],
			    minHeight: 200,
			    height: 300,
			    bodyPadding: 0,
			    split: true,
			    items:[clusterPanel]
			},{ region: 'center',
			    layout: 'fit',
			    id: 'expressionDetailsPane',
			    minHeight: 100,
			    collapsible: false,
//			    headerPosition: 'bottom',
			    title: 'Expression Details',
			    tools: [
				{ type:'search',
				  tooltip: 'Search for genes matching the current consensus pattern',
				  handler: function(e, el,o,t) {
				      detailPanel.searchSimilar();
				  }
				},{ type:'collapse',
				  tooltip: 'Run GO enrichment analysis on the current gene set',
				  handler: function(e, el,o,t) {
				      if(detailPanel.hasOwnProperty('genecldata')) {
					  // write out current gene set
					  geinfostore.getProxy().setExtraParam("genes",JSON.stringify(detailPanel.genecldata.matrix.rows))
					  geinfostore.load();
					  infotab.child("#enrichmenttab").tab.show();
					  // show the info tab
					  infotab.setActiveTab(3);
				      }
				  }
				},{ type:'gear',
				  tooltip: 'Settings',
				  handler: function(e, el,o,t) {
				      detailPanelGearMenu.showBy(t);
				  }
				},
				{ type:'save',
 				  tooltip: 'Save image as SVG file',
				  handler: function(e,el,o,t) {
				      var svg=d3.select("#geneclsvg");
				      if(!svg.empty()) {
					  svg.selectAll(".geneRowLabel").style("dominant-baseline","middle").style("font","12px sans-serif");;
					  //var b64 = Base64.encode(el.html());
					  var b64 = window.btoa(svg.node().parentNode.innerHTML);
					  writeAndClickLink('data:application/octet-stream;base64,\n'+b64,'genes.svg')
				      }
				  }
				}
			    ],
			    header: true,
			    items:[detailPanel],
			    autoScroll: true,
			    autoShow: true,
			    /*listeners: {
				afterrender: function(panel) {
				    console.log("boo");
				    var header=panel.getHeader();
				    header.insert(1,[ngenesSlider]);
				}
			    }*/
			}]
		},{ region: 'east',
		    collapsible: true,
		    title: 'Info',
		    split: true,
		    layout: 'fit',
		    width: '30%',
		    minWidth: 100,
		    minHeight: 140,
		    bodyPadding: 0,
		    items:[infotab]
		}]
	});


    if(Ext.util.Cookies.get("hidetutorial")==null) {
	tutorialWindow.show(); 
    }


});



// quick helper function to provide an internal download link
function writeAndClickLink(url,download) {
    var link=d3.select("body").append('a');
    link.attr('href',url).attr('download',download);
    var evObj = document.createEvent('MouseEvents');
    evObj.initMouseEvent('click', true, true);
    link.node().dispatchEvent(evObj);
    $(link.node()).remove();
}
