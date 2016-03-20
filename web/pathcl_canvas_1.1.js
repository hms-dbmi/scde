Ext.require(['*']);

d3.selection.prototype.moveToFront = function() {
  return this.each(function(){
  this.parentNode.appendChild(this);
  });
};

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

    
/* 2D embedding panel */
    var embeddingPanel = Ext.create('Ext.panel.Panel', {
	layout: 'fit',
	bodyPadding:0,
	autoScroll: false,
	draw: function(cols) {
	    if(!clusterPanel.pathcldata.hasOwnProperty('embedding')) {
		return;
	    }
	    if(cols===undefined) {
		// use last colcol row for colors
		cols=[];
		if(clusterPanel.pathcldata.hasOwnProperty('colcols')) {
		    var i=clusterPanel.pathcldata.colcols.dim[0]-1;
		    for(var j=0; j<clusterPanel.pathcldata.colcols.dim[1]; j++) { 
			cols.push(clusterPanel.pathcldata.colcols.data[i*clusterPanel.pathcldata.colcols.dim[1]+j])
		    }
		} else {
		    // color everything gray
		    for(var j=0;j<clusterPanel.pathcldata.embedding.dim[0];j++) {
			cols.push("#808080")
		    }
		}
	    }
	    
	    embeddingPanel.cols=cols;
	    embeddingPanel.redraw();
	},
	recolor: function(cols) {
	    if(embeddingPanel.hasOwnProperty('cols')) {
		// recolor
		embeddingPanel.cols=cols;
		d3.select("#embedding").selectAll("circle").attr("fill",function(d,i) {
		    return(embeddingPanel.cols[d3.select(this).attr("celli")-1])
		})
	    } else { // just redraw everything
		embeddingPanel.draw(cols);
	    }
	},
	redraw: function() {
	    if(!clusterPanel.hasOwnProperty('pathcldata') || !clusterPanel.pathcldata.hasOwnProperty('embedding') || !embeddingPanel.hasOwnProperty('cols')) {
		return;
	    }

	    var padding=12;
	    var el = d3.select(embeddingPanel.getLayout().getElementTarget().dom)
	    var s=embeddingPanel.getSize();
	    d3.select("#embedding").remove();
	    var svg = el.append("svg").attr("id","embedding").attr("width",s.width+"px").attr("height",s.height+"px").attr('xmlns','http://www.w3.org/2000/svg');

	    var xScale=d3.scale.linear().domain(clusterPanel.pathcldata.embedding.xrange).range([padding,s.width-2*padding]);
	    var yScale=d3.scale.linear().domain(clusterPanel.pathcldata.embedding.yrange).range([padding,s.height-2*padding]);

	    //var yAxis=d3.svg.axis().scale(yScale).orient("left");
	    //var xAxis=d3.svg.axis().scale(xScale).orient("bottom");
	    //svg.append("g").attr("class","axis").attr("transform", "translate(0," + (s.height - 2*padding) + ")").call(xAxis);
	    //svg.append("g").attr("class","axis").attr("transform", "translate(" + 2*padding + ",0)").call(yAxis);
	    //svg.append("g").attr("class","axis").call(xAxis);
	    embeddingPanel.embarray = $.map(clusterPanel.pathcldata.embedding.data, function(value, index) {
		value.push(index); return [value];
	    });
	    svg.selectAll("circle")
		.data(embeddingPanel.embarray)
		.enter()
		.append("circle")
		.attr("id",function(d) {return("cell"+(d[2]-1))})
		.attr("celli",function(d) {return(d[2])})
		.attr("cx",function(d) { return(xScale(d[0])) })
		.attr("cy",function(d) { return(yScale(d[1])) })
		.attr("fill",function(d) {return(embeddingPanel.cols[d[2]-1])})
		.attr("r",5).on("mouseover",function() {
		    d3.selectAll("#embedding circle.selected").classed("selected",false);
		    var rx=d3.select(this).classed("selected",true).attr("celli")-1;
		    var v=clusterPanel.s;
		    clusterPanel.evctx.clearRect(0,0,v.width,v.height);
		    clusterPanel.evctx.beginPath();
		    var mx=hc.hmleft()+(rx+0.5)*hc.hmwidth(v.width)/clusterPanel.pathcldata.matrix.dim[1];
		    clusterPanel.evctx.moveTo(mx,hc.cctop(clusterPanel.pathcldata))
		    clusterPanel.evctx.lineTo(mx,v.height-hc.margins.bottom)
		    clusterPanel.evctx.stroke();
		    if(detailPanel.hasOwnProperty('evctx')) {
			v=detailPanel.getSize();
			detailPanel.evctx.clearRect(0,0,v.width,v.height);
			detailPanel.evctx.beginPath();
			detailPanel.evctx.moveTo(mx,hc.cctop(detailPanel.genecldata))
			detailPanel.evctx.lineTo(mx,v.height-hc.margins.bottom)
			detailPanel.evctx.stroke();
		    }
		    
		}).append("svg:title").text(function(d){ return(d[3]) });
	    
	    svg.on('mouseout', function() {
		var v=clusterPanel.s;
		clusterPanel.evctx.clearRect(0,0,v.width,v.height);
		v=detailPanel.getSize();
		if(detailPanel.hasOwnProperty('evctx')) {
		    detailPanel.evctx.clearRect(0,0,v.width,v.height);
		}
	    });
		
	},
	listeners: {
	    resize: function(cmp,width,height,oldWidth,oldHeight,opts) {
		cmp.redraw();
	    }
	}
    });
    

    /* DETAILED CLUSTERING VIEW */
    var detailMode=1; // 1: gene 2: pathway
    var detailItemList={}; // pathway or gene list
    var detailNGenes=20; // max genes to show 

    var detailPanel = Ext.create('Ext.panel.Panel', {
	layout: 'fit',
	bodyPadding: 5,
	autoScroll: false,
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
		    if(data.matrix.dim[0]==1) data.matrix.rows=[data.matrix.rows]
		    detailPanel.genecldata=data;
		    detailPanel.mode=1;
		    Ext.getCmp('expressionDetailsPane').setTitle("Expression Details: selected genes");
		    detailPanel.redraw(detailPanel.genecldata)
		    detailPanelGearMenu.getComponent(0).disable()
		    detailPanelGearMenu.getComponent(1).enable()
		    detailPanelGearMenu.getComponent(3).disable()
		    detailPanel.setLoading(false);
		}
	    });
	},
	showPathways: function(ids) {
	    if(ids === undefined || ids.length==0) return;
	    detailPanel.setLoading(true);
	    var ngenes=detailPanelGearMenu.getComponent(0).getValue();
	    var twosided=detailPanelGearMenu.getComponent(3).checked;
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
		    if(data.matrix.dim[0]==1) data.matrix.rows=[data.matrix.rows]
		    detailPanel.genecldata=data;
		    detailPanel.mode=2; 
		    Ext.getCmp('expressionDetailsPane').setTitle("Expression Details: top genes in specified pathways");
		    detailPanel.redraw(detailPanel.genecldata)
		    detailPanelGearMenu.getComponent(0).enable()
		    detailPanelGearMenu.getComponent(1).enable()
		    detailPanelGearMenu.getComponent(3).enable()
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
	    var twosided=detailPanelGearMenu.getComponent(3).checked;
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
		    if(data.matrix.dim[0]==1) data.matrix.rows=[data.matrix.rows]
		    detailPanel.genecldata=data;
		    detailPanel.mode=3;
		    Ext.getCmp('expressionDetailsPane').setTitle("Expression Details: genes matching specified pattern");
		    detailPanel.redraw(detailPanel.genecldata)
		    detailPanelGearMenu.getComponent(0).enable()
		    detailPanelGearMenu.getComponent(1).enable()
		    detailPanelGearMenu.getComponent(3).enable()
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
	    
	    
	    $('#genecl').remove();
	    delete detailPanel.evctx;
	    $('#geneclev').remove();
	    $('#gclCD').remove();
	    var s=detailPanel.getSize();
	    // adjust height to match maxRowHeight (=15) if needed
	    if((s.height-hc.hmtop(data)-hc.margins.bottom)/data.matrix.dim[0] > 15) {
		s.height=15*data.matrix.dim[0]+hc.margins.bottom+hc.hmtop(data);
	    }
	    detailPanel.s=s;
	    
	    detailPanel.body.update('<div id="gclCD" style="padding:0 0 0 0; position:relative"><canvas id="genecl" width='+s.width+' height='+s.height+' style="background-color:transparent; position:absolute; left: 0; top: 0; z-index: 0;"></canvas> <canvas id="geneclev" width='+s.width+' height='+s.height+' style="background-color:transparent; position:absolute; left: 0; top: 0; z-index: 1;"></canvas></div>')
	    
	    var gctx= $('#genecl')[0].getContext('2d');
	    if(data.hasOwnProperty('colcols')) {
		//colcols
		drawHeatmap(gctx,data.colcols,hc.hmleft(),hc.cctop(data),hc.hmwidth(s.width), hc.ccheight(data));
		//rowcols
		drawHeatmap(gctx,data.rowcols,hc.margins.left,hc.hmtop(data),hc.rowcolWidth,hc.hmheight(data,s.height),false)
	    }

	    //main heatmap
	    drawHeatmap(gctx,data.matrix,hc.hmleft(),hc.hmtop(data),hc.hmwidth(s.width),hc.hmheight(data,s.height),true)
	    // event handling
	    var evc=document.getElementById("geneclev");
	    var evctx= evc.getContext('2d');
	    detailPanel.evctx=evctx;
	    try {
		evctx.setLineDash([5]);
	    } catch (err) {}
	    evctx.fillStyle='black'; 
	    evctx.font="bold 15px Arial";
	    
	    var evtRect = evc.getBoundingClientRect();


	    evc.addEventListener('click', function(evt) {
		if(!clusterPanel.hasOwnProperty("pathcldata") || !clusterPanel.pathcldata.hasOwnProperty('embedding')) { return; }
		var mx=evt.clientX-evtRect.left; var my=evt.clientY-evtRect.top;
		if(my<=s.height-hc.margins.bottom && mx<=s.width-hc.margins.right && mx>=hc.hmleft() && my>=hc.cctop(data)) { 
		    if(my>=hc.hmtop(data)) {
			var ry=Math.floor((my-hc.hmtop(data))/hc.hmheight(data,s.height)*data.matrix.dim[0]);
			// update embedding
			embeddingPanel.recolor(heatmapRowColors(data.matrix,ry));
			Ext.getCmp("embeddingDiv").setTitle("2D Embedding: gene "+data.matrix.rows[ry])
		    } else {
			var ry=Math.floor((my-hc.cctop(data))/hc.ccheight(data)*data.colcols.dim[0]);
			embeddingPanel.recolor(heatmapRowColors(data.colcols,ry));
			Ext.getCmp("embeddingDiv").setTitle("2D Embedding: consensus pattern")
		    }
		}
	    }, false);
	    
	    evc.addEventListener('mousemove', function(evt) {
		var mx=evt.clientX-evtRect.left; var my=evt.clientY-evtRect.top;

		evctx.clearRect(0,0,s.width,s.height);
		var v=clusterPanel.s;
		clusterPanel.evctx.clearRect(0,0,v.width,v.height);
		if(my<=s.height-hc.margins.bottom && mx<=s.width-hc.margins.right && mx>=hc.hmleft() && my>=hc.cctop(data)) { 
		    //mx=hc.hmleft()+(rx+0.5)*hc.hmwidth(s.width)/data.matrix.dim[1];
		    evctx.beginPath(); 
		    evctx.moveTo(mx,hc.cctop(data)); evctx.lineTo(mx,s.height-hc.margins.bottom);

		    var rx=Math.floor((mx-hc.hmleft())/hc.hmwidth(s.width)*data.matrix.dim[1]);
		    
		    // update the line in the cluster panel
		    clusterPanel.evctx.beginPath();
		    clusterPanel.evctx.moveTo(mx,hc.cctop(clusterPanel.pathcldata))
		    clusterPanel.evctx.lineTo(mx,v.height-hc.margins.bottom)
		    clusterPanel.evctx.stroke();

		    if(mx>s.width/2) { 
			evctx.textAlign="end"; mx-=5; 
		    } else {
			evctx.textAlign="start"; mx+=5;
		    }
		    evctx.textBaseline="bottom";
		    evctx.fillText("cell: "+data.matrix.cols[rx],mx,my-3); 
		    


		    d3.selectAll("#embedding circle.selected").classed("selected",false);
		    d3.select("#cell"+rx).classed("selected",true).moveToFront();

		    if(my>=hc.hmtop(data)) {
			var ry=Math.floor((my-hc.hmtop(data))/hc.hmheight(data,s.height)*data.matrix.dim[0]);
			//my=hc.hmtop(data)+(ry+0.5)*hc.hmheight(data,s.height)/data.matrix.dim[0];
			evctx.moveTo(hc.margins.left,my); evctx.lineTo(s.width-hc.margins.right,my);

			evctx.textBaseline="top";
			evctx.fillText("gene: "+data.matrix.rows[ry],mx,my+3); 	
		    } else {
			evctx.moveTo(hc.hmleft(),my); evctx.lineTo(s.width-hc.margins.right,my);
		    }
		    evctx.stroke();

		}
	    }, false);

	    evc.addEventListener('mouseout', function(evt) {
		evctx.clearRect(0,0,s.width,s.height);
		var v=clusterPanel.s;
		clusterPanel.evctx.clearRect(0,0,v.width,v.height);
		d3.selectAll("#embedding circle.selected").classed("selected",false);
	    });

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
            {text: "Variance", width: 100, dataIndex: 'var', sortable: true}
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
			change: {buffer: 600, fn: function(field, value) {
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
			change: {buffer: 400, fn: function(field, value) {
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

    
    // draw dendrogram in a given context
    // ctx - canvas 2d context
    // d - dendrogram structure (t(merge), order, height)
    // x,y,width,height - placement
    function drawDendrogram(ctx, d, x, y, width, height) {
	var nmerges=d.merge.length;
	var xstep=width/d.order.length;
	var lp,lpy,rp,rpy,jx,jy;
	var maxHeight=Math.max.apply(null, d.height);
	var heightScale=height/maxHeight;
	ctx.beginPath();

	for(var i=0, mergex=[]; i<nmerges; i++) {
	    jy=y+(maxHeight-d.height[i])*heightScale;
	    lp=d.merge[i]; // left position (starts as an index)
	    if(lp<0) { // leaf
		lp=d.order[-1*lp-1]; 
		lp=x+(lp-0.5)*xstep;
		lpy=y+height;
	    } else { // inner node
		lpy=y+(maxHeight-d.height[lp-1])*heightScale; lp=mergex[lp-1]; 
	    }
	    i++; rp=d.merge[i]; // right position
	    if(rp<0) { // leaf
		rp=d.order[-1*rp-1]; 
		rp=x+(rp-0.5)*xstep;
		rpy=y+height;
	    } else { // inner node
		rpy=y+(maxHeight-d.height[rp-1])*heightScale; rp=mergex[rp-1]; 
	    }
	    jx=(lp+rp)/2;  jy=y+(maxHeight-d.height[(i-1)/2])*heightScale;
	    mergex.push(jx);
	    ctx.moveTo(lp,lpy); ctx.lineTo(lp,jy); ctx.lineTo(rp,jy); ctx.lineTo(rp,rpy);
	    
	}
	ctx.lineWidth=1.5; ctx.lineJoin='miter'; ctx.strokeStyle='black'; ctx.stroke();
    }

    // get heatmap row (in colors)
    function heatmapRowColors(d,i) {
	var cols=[];
	if(d.hasOwnProperty('colors')) {
	    // establish mapping
	    if(d.hasOwnProperty('zlim')) {
		zlim=d.zlim;
	    } else {
		zlim.push(-1*Math.max.apply(null,d.data.map(Math.abs)));
		zlim.push(-1*zlim[0]);
	    }
	    zlim.push(d.colors.length/(zlim[1]-zlim[0])); // zlim step
	    var val;
	    for(var j=0; j<d.dim[1];j++) { // columns 
		val=(d.data[i*d.dim[1]+j] - zlim[0])*zlim[2];
		if(val<0) { val=0; } else if(val>=d.colors.length) { val=d.colors.length-1; }
		val=d.colors[Math.floor(val)];
		cols.push(val)
	    }
	} else {
	    // interpret the values as colors
	    for(var j=0; j<d.dim[1];j++) { // columns 
		cols.push(d.data[i*d.dim[1]+j])
	    }
	}
	return(cols);
    }

     // draw heatmap
    function drawHeatmap(ctx,d,x,y,width,height,rowNames,maxRowHeight,maxFontSize,minFontSize) {
	rowNames = typeof rowNames !== 'undefined' ? rowNames : false;
	maxFontSize = typeof maxFontSize !== 'undefined' ? maxFontSize : 100;
	minFontSize = typeof minFontSize !== 'undefined' ? minFontSize : 6;
	maxRowHeight = typeof maxRowHeight !== 'undefined' ? maxRowHeight : 100;
	if(height>d.dim[0]*maxRowHeight) { height=d.dim[0]*maxRowHeight;}
	var rowHeight=height/d.dim[0]; var colWidth=width/d.dim[1];
	var mC=d.hasOwnProperty('colors'); // perform color mapping on the fly
	if(rowNames && !d.hasOwnProperty('rows')) { rowNames=false; }
	if(rowNames) { 
	    var fontSize=(rowHeight*0.95); 
	    if(fontSize>maxFontSize) { fontSize=maxFontSize; };
	    if(fontSize>=minFontSize) { 
		ctx.font=fontSize+"px Arial"; ctx.textAlign='left'; ctx.textBaseline="middle";
	    } else { 
		rowNames=false;
	    }
	}
	var zlim=[];
	if(mC) {
	    if(d.hasOwnProperty('zlim')) {
		zlim=d.zlim;
	    } else {
		zlim.push(-1*Math.max.apply(null,d.data.map(Math.abs)));
		zlim.push(-1*zlim[0]);
	    }
	    zlim.push(d.colors.length/(zlim[1]-zlim[0])); // zlim step
	}
	var val=0;
	for(var i=0; i<d.dim[0]; i++) { // rows
	    for(var j=0; j<d.dim[1]; j++) { // columns
		if(mC) { // colors were specified
		    val=(d.data[i*d.dim[1]+j] - zlim[0])*zlim[2];
		    if(val<0) { val=0; } else if(val>=d.colors.length) { val=d.colors.length-1; }
		    val=d.colors[Math.floor(val)];
		} else { // interpret values as colors directly
		    val=d.data[i*d.dim[1]+j]
		}
		ctx.fillStyle=val; ctx.fillRect(x+j*colWidth,y+i*rowHeight,colWidth,rowHeight);

	    }
	    if(rowNames) { ctx.fillStyle='black'; ctx.fillText(d.rows[i],x+width+3,y+(i+0.5)*rowHeight); }
	}
	ctx.lineWidth=1; ctx.strokeStyle='black'; ctx.strokeRect(x,y,width,height);
	//if(zlim!==undefined) {return(zlim[1])};
    }

    function updatePathclInfo(pathcl) {
	if(pathcl == currentPathCl) return;
	clinfostore.getProxy().setExtraParam("pathcl",pathcl)
	clinfostore.load();
	infotab.child("#clustertab").tab.show();
	infotab.setActiveTab(2);
	infotab.getActiveTab().setTitle("Aspect "+pathcl)
	currentPathCl=pathcl;
    }
    
    /* heatmap config */
    var hc = {
	spacing: 5,
	dendSpacing: 1,
	geneUnitHeight: 15,
	margins: {top:2,right:80,bottom:15,left:1},
	colcolUnitHeight: 10,
	rowcolWidth: 10,
	colDendHeight: 50,
	trim: 0,
	// heatmap left position
	hmleft: function() { return(this.margins.left+this.rowcolWidth+this.spacing)},
	// colcol top position
	cctop: function(data) { return(data.hasOwnProperty('coldend')? this.margins.top+this.colDendHeight+this.dendSpacing : this.margins.top) },
	// colcol height
	ccheight: function(data) { return(this.colcolUnitHeight*data.colcols.dim[0]) },
	// heatmap top position
	hmtop: function(data) { 
	    if(data.hasOwnProperty('colcols')) {
		return(data.hasOwnProperty('coldend')? this.margins.top+this.colDendHeight+this.spacing+this.dendSpacing + this.ccheight(data) : this.margins.top+this.spacing+this.dendSpacing + this.ccheight(data))
	    } else {
		return(data.hasOwnProperty('coldend')? this.margins.top+this.colDendHeight+this.spacing : this.margins.top)
	    }
	},
	// heatmap hight
	hmheight: function(data,totalHeight) {
	    return((totalHeight-this.margins.bottom) - this.hmtop(data));
	},
	// heatmap width
	hmwidth: function(totalWidth) { return(totalWidth-this.margins.left-this.rowcolWidth-this.spacing-this.margins.right)},
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
	    Ext.Ajax.request({
		url: 'pathcl.json',
		method: 'POST',          
		waitTitle: 'Connecting',
		waitMsg: 'Sending data...',                                     
		scope:this,
		failure: function(r,o){console.log('failure:'); console.log(r);},
		success: function(response) {
		    var data = Ext.JSON.decode(response.responseText)
		    if(data.matrix.dim[0]==1) data.matrix.rows=[data.matrix.rows]
		    if(data.colcols.dim[0]==1) data.colcols.rows=[data.colcols.rows]
		    clusterPanel.pathcldata=data;
		    clusterPanelGearMenu.getComponent(0).suspendEvents();
		    clusterPanelGearMenu.getComponent(0).setValue(data.matrix.zlim[1]);
		    //clusterPanelGearMenu.getComponent(0).setMaxValue(data.matrix.range[1]);
		    clusterPanelGearMenu.getComponent(0).resumeEvents()
		    if(data.hasOwnProperty('trim')) {
			hc.trim=data.trim;
			detailPanelGearMenu.getComponent(1).suspendEvents();
			detailPanelGearMenu.getComponent(1).setValue(data.trim);
			detailPanelGearMenu.getComponent(1).resumeEvents();
		    }
		    clusterPanel.redraw(clusterPanel.pathcldata);
		    if(data.hasOwnProperty('embedding')) {
			Ext.getCmp("embeddingDiv").show();
			embeddingPanel.draw();
		    }
		    clusterPanel.setLoading(false);
		}
	    })
	},
	redraw: function(data) {
	    if(data === undefined) {
		if(clusterPanel.hasOwnProperty('pathcldata')) { data=clusterPanel.pathcldata; } else { clusterPanel.reload(); return; }
	    }
	    
	    $('#pathcl').remove();
	    $('#pathclev').remove();
	    delete clusterPanel.evctx;
	    $('#pclCD').remove();
	    //clusterPanel.update("");
	    var s=clusterPanel.getSize();
	    clusterPanel.s=s;
	    clusterPanel.body.update('<div id="pclCD" style="padding:0 0 0 0; position:relative"><canvas id="pathcl" width='+s.width+' height='+s.height+' style="background-color:transparent; position:absolute; left: 0; top: 0; z-index: 0;"></canvas> <canvas id="pathclev" width='+s.width+' height='+s.height+' style="background-color:transparent; position:absolute; left: 0; top: 0; z-index: 1;"></canvas></div>')
	    
	    
            var ctx= $('#pathcl')[0].getContext('2d');
	    drawDendrogram(ctx,data.coldend,hc.hmleft(),hc.margins.top,hc.hmwidth(s.width),hc.colDendHeight);
	    
	    //colcols
	    drawHeatmap(ctx,data.colcols,hc.hmleft(),hc.cctop(data),hc.hmwidth(s.width), hc.ccheight(data));
	    
	    //rowcols
	    drawHeatmap(ctx,data.rowcols,hc.margins.left,hc.hmtop(data),hc.rowcolWidth,hc.hmheight(data,s.height))

	    //main heatmap
	    drawHeatmap(ctx,data.matrix,hc.hmleft(),hc.hmtop(data),hc.hmwidth(s.width),hc.hmheight(data,s.height),rowNames=true)
	    
	    // event handling
	    var evc=document.getElementById("pathclev");
	    var evctx= evc.getContext('2d');
	    clusterPanel.evctx=evctx;

	    //evctx.strokeRect(0,0,s.width-10,s.height-10);
	    try {
		evctx.setLineDash([5]);
	    } catch (err) {}
	    evctx.fillStyle='black'; 
	    evctx.font="bold 15px Arial";

	    var evtRect = evc.getBoundingClientRect();
	    
	    evc.addEventListener('click', function(evt) {
		var mx=evt.clientX-evtRect.left; var my=evt.clientY-evtRect.top;
		if(my<=s.height-hc.margins.bottom && mx<=s.width-hc.margins.right && mx>=hc.hmleft() && my>=hc.cctop(data)) { 
		//if(my<=s.height-hc.margins.bottom && mx<=s.width-hc.margins.right && mx>=hc.hmleft()) { 
		    if(my>=hc.hmtop(data)) {
			var ry=Math.floor((my-hc.hmtop(data))/hc.hmheight(data,s.height)*data.matrix.dim[0]);
			evctx.strokeStyle='red'; evctx.lineWidth=2;
			try { evctx.setLineDash([0]); } catch(err) {}
			//evctx.strokeRect(hc.hmleft(),hc.hmtop(data)+ry*hc.hmheight(data,s.height)/data.matrix.dim[0],hc.hmwidth(s.width),hc.hmheight(data,s.height)/data.matrix.dim[0]); 
			evctx.strokeRect(0,hc.hmtop(data)+ry*hc.hmheight(data,s.height)/data.matrix.dim[0],s.width,hc.hmheight(data,s.height)/data.matrix.dim[0]); 
			evctx.strokeStyle='black'; evctx.lineWidth=1;
			try { evctx.setLineDash([5]); } catch(err) {}
			updatePathclInfo(data.matrix.rows[ry])

			// update embedding
			
			if(data.hasOwnProperty('embedding')) {
			    embeddingPanel.recolor(heatmapRowColors(data.matrix,ry));
			    Ext.getCmp("embeddingDiv").setTitle("2D Embedding: aspect "+data.matrix.rows[ry]);
			}
		    } else {
			if(data.hasOwnProperty('embedding')) {
			    var ry=Math.floor((my-hc.cctop(data))/hc.ccheight(data)*data.colcols.dim[0]);
			    embeddingPanel.recolor(heatmapRowColors(data.colcols,ry));
			    Ext.getCmp("embeddingDiv").setTitle("2D Embedding: metadata "+data.colcols.rows[ry])
			}
			
		    }
		}
	    }, false);


	    evc.addEventListener('mousemove', function(evt) {
		var mx=evt.clientX-evtRect.left; var my=evt.clientY-evtRect.top;

		evctx.clearRect(0,0,s.width,s.height);
		if(detailPanel.hasOwnProperty('evctx')) {
		    var v=detailPanel.s;
		    detailPanel.evctx.clearRect(0,0,v.width,v.height);
		}
		if(my<=s.height-hc.margins.bottom && mx<=s.width-hc.margins.right && mx>=hc.hmleft() && my>=hc.cctop(data)) { 
		    //mx=hc.hmleft()+(rx+0.5)*hc.hmwidth(s.width)/data.matrix.dim[1];
		    evctx.beginPath(); 
		    evctx.moveTo(mx,hc.cctop(data)); evctx.lineTo(mx,s.height-hc.margins.bottom);
		    // update the line in the gene panel
		    if(detailPanel.hasOwnProperty('evctx')) {
			detailPanel.evctx.beginPath();
			detailPanel.evctx.moveTo(mx,hc.cctop(detailPanel.genecldata))
			detailPanel.evctx.lineTo(mx,v.height-hc.margins.bottom)
			detailPanel.evctx.stroke();
		    }
		    var rx=Math.floor((mx-hc.hmleft())/hc.hmwidth(s.width)*data.matrix.dim[1]);
		    if(mx>s.width/2) { 
			evctx.textAlign="end"; mx-=5; 
		    } else {
			evctx.textAlign="start"; mx+=5;
		    }
		    evctx.textBaseline="bottom";
		    evctx.fillText("cell: "+data.matrix.cols[rx],mx,my-3); 

		    if(data.hasOwnProperty('embedding')) {
			d3.selectAll("#embedding circle.selected").classed("selected",false);
			d3.select("#cell"+rx).classed("selected",true).moveToFront();
		    }
		    
		    if(my>=hc.hmtop(data)) {
			var ry=Math.floor((my-hc.hmtop(data))/hc.hmheight(data,s.height)*data.matrix.dim[0]);
			//my=hc.hmtop(data)+(ry+0.5)*hc.hmheight(data,s.height)/data.matrix.dim[0];
			evctx.moveTo(hc.margins.left,my); evctx.lineTo(s.width-hc.margins.right,my);
			
			//evctx.fillText("cell: "+data.matrix.cols[rx],mx+5,my-5); 
			evctx.textBaseline="top";
			evctx.fillText("aspect: "+data.matrix.rows[ry],mx,my+3);

		    } else {
			evctx.moveTo(hc.hmleft(),my); evctx.lineTo(s.width-hc.margins.right,my);
			var ry=Math.floor((my-hc.cctop(data))/hc.ccheight(data)*data.colcols.dim[0]);
			var val=data.colcols.rows[ry]; 
			if(val!==undefined) {
			    evctx.textBaseline="top";
			    evctx.fillText("metadata: "+val,mx,my+3);
			}
		    }
		    evctx.stroke();
		}
	    }, false);
	    

	    evc.addEventListener('mouseout', function(evt) {
		evctx.clearRect(0,0,s.width,s.height);
		var v=detailPanel.getSize();
		if(detailPanel.hasOwnProperty('evctx')) {
		    detailPanel.evctx.clearRect(0,0,v.width,v.height);
		}
		d3.selectAll("#embedding circle.selected").classed("selected",false);
	    });
	    
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
			clusterPanel.pathcldata.matrix.zlim=[-1*v,v];
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
			    title: 'Pathway Overdispersion',
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
				  tooltip: 'Save image',
				  handler: function(e,el,o,t) {
				      if($('#pathcl').length==0) { return; }
				      var changingImage = Ext.create('Ext.Img', {
					  src: $('#pathcl')[0].toDataURL("image/png",1.0)
				      });
				      win = new Ext.Window({
					  title: 'exported image: use right click to save the image',
					  layout: 'fit',
					  autoScroll: true,
					  modal: true,
					  closeAction: 'hide',
					  items:[changingImage]
				      });
				      win.show();
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
				  tooltip: 'Save image',
				  handler: function(e,el,o,t) {
				      if($('#genecl').length==0) { return; }
				      var changingImage = Ext.create('Ext.Img', {
					  src: $('#genecl')[0].toDataURL("image/png",1.0)
				      });
				      win = new Ext.Window({
					  title: 'exported image: use right click to save the image',
					  layout: 'fit',
					  autoScroll: true,
					  modal: true,
					  closeAction: 'hide',
					  items:[changingImage]
				      });
				      win.show();
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
		    split: true,
		    layout: 'border',
		    width: '30%',
		    title:'Info',
		    bodyPadding: 0,
		    items:[{
			region:'center',
			layout:'fit',
			minWidth: 100,
			minHeight: 140,
			bodyPadding: 0,
			items:[infotab]
		    },{
			title:'2D Embedding',
			id:'embeddingDiv',
			region:'south',
			layout: 'fit',
			minWidth: 100,
			minHeight: 100,
			height:'40%',
			hidden:true,
			collapsible: true,
			split:true,
			bodyPadding: 0,
			tools:[{ type:'save',
				 tooltip: 'Save image',
				 handler: function(e,el,o,t) {
				     var svge=d3.select("#embedding");
				     if(!svge.empty()) {
					 var html=svge.attr("version", 1.1)
					     .attr("xmlns", "http://www.w3.org/2000/svg")
					     .node().parentNode.innerHTML;
					 var imgsrc = 'data:image/svg+xml;base64,'+ btoa(html);
					 win = new Ext.Window({
					     title: 'exported image: use the link below save the image',
					     layout: 'fit',
					     html: '<a href='+imgsrc+'>link</a>'
					 });
					 win.show();
					 
				     }
				 }
			       }],
			items:[embeddingPanel]
		    }]
		}]
	});


    if(Ext.util.Cookies.get("hidetutorial")==null) {
	tutorialWindow.show(); 
    }


});
// google analytics code
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-33018606-2', 'auto');
  ga('send', 'pageview');

