(function(e){function t(t){for(var s,n,a=t[0],l=t[1],c=t[2],d=0,v=[];d<a.length;d++)n=a[d],Object.prototype.hasOwnProperty.call(o,n)&&o[n]&&v.push(o[n][0]),o[n]=0;for(s in l)Object.prototype.hasOwnProperty.call(l,s)&&(e[s]=l[s]);u&&u(t);while(v.length)v.shift()();return i.push.apply(i,c||[]),r()}function r(){for(var e,t=0;t<i.length;t++){for(var r=i[t],s=!0,a=1;a<r.length;a++){var l=r[a];0!==o[l]&&(s=!1)}s&&(i.splice(t--,1),e=n(n.s=r[0]))}return e}var s={},o={app:0},i=[];function n(t){if(s[t])return s[t].exports;var r=s[t]={i:t,l:!1,exports:{}};return e[t].call(r.exports,r,r.exports,n),r.l=!0,r.exports}n.m=e,n.c=s,n.d=function(e,t,r){n.o(e,t)||Object.defineProperty(e,t,{enumerable:!0,get:r})},n.r=function(e){"undefined"!==typeof Symbol&&Symbol.toStringTag&&Object.defineProperty(e,Symbol.toStringTag,{value:"Module"}),Object.defineProperty(e,"__esModule",{value:!0})},n.t=function(e,t){if(1&t&&(e=n(e)),8&t)return e;if(4&t&&"object"===typeof e&&e&&e.__esModule)return e;var r=Object.create(null);if(n.r(r),Object.defineProperty(r,"default",{enumerable:!0,value:e}),2&t&&"string"!=typeof e)for(var s in e)n.d(r,s,function(t){return e[t]}.bind(null,s));return r},n.n=function(e){var t=e&&e.__esModule?function(){return e["default"]}:function(){return e};return n.d(t,"a",t),t},n.o=function(e,t){return Object.prototype.hasOwnProperty.call(e,t)},n.p="";var a=window["webpackJsonp"]=window["webpackJsonp"]||[],l=a.push.bind(a);a.push=t,a=a.slice();for(var c=0;c<a.length;c++)t(a[c]);var u=l;i.push([0,"chunk-vendors"]),r()})({0:function(e,t,r){e.exports=r("56d7")},"29bc":function(e,t,r){},"2f3f":function(e,t,r){},"56d7":function(e,t,r){"use strict";r.r(t);r("e260"),r("e6cf"),r("cca6"),r("a79d");var s=r("2b0e"),o=(r("d3b7"),r("bc3a")),i=r.n(o),n={},a=i.a.create(n);a.interceptors.request.use((function(e){return e}),(function(e){return Promise.reject(e)})),a.interceptors.response.use((function(e){return e}),(function(e){return Promise.reject(e)})),Plugin.install=function(e){e.axios=a,window.axios=a,Object.defineProperties(e.prototype,{axios:{get:function(){return a}},$axios:{get:function(){return a}}})},s["a"].use(Plugin);Plugin;var l=function(){var e=this,t=e.$createElement,r=e._self._c||t;return r("v-app",[r("nav",[r("v-app-bar",{attrs:{color:"primary",dark:""}},[r("v-app-bar-nav-icon",{on:{click:function(t){e.drawer=!0}}}),r("v-toolbar-title",[e._v("LUMC Mutalyzer 3")])],1),r("v-navigation-drawer",{attrs:{absolute:"",temporary:""},model:{value:e.drawer,callback:function(t){e.drawer=t},expression:"drawer"}},[r("v-list",{attrs:{nav:"",dense:""}},[r("v-list-item-group",{attrs:{"active-class":"primary--text text--accent-4"}},[r("v-list-item",{attrs:{to:"/"}},[r("v-list-item-icon",[r("v-icon",[e._v("mdi-home")])],1),r("v-list-item-title",[e._v("Home")])],1),r("v-list-item",{attrs:{to:"/namechecker"}},[r("v-list-item-icon",[r("v-icon",[e._v("mdi-description")])],1),r("v-list-item-title",[e._v("Name Checker")])],1),r("v-list-item",{attrs:{to:"/positionconverter"}},[r("v-list-item-icon",[r("v-icon",[e._v("mdi-description")])],1),r("v-list-item-title",[e._v("Position Converter")])],1),r("v-list-item",{attrs:{to:"/descriptionextractor"}},[r("v-list-item-icon",[r("v-icon",[e._v("mdi-description")])],1),r("v-list-item-title",[e._v("Description Extractor")])],1),r("v-list-item",{attrs:{to:"/about"}},[r("v-list-item-icon",[r("v-icon",[e._v("mdi-description")])],1),r("v-list-item-title",[e._v("About")])],1)],1)],1)],1)],1),r("v-main",[r("router-view")],1)],1)},c=[],u={name:"App",data:function(){return{drawer:!1}}},d=u,v=r("2877"),p=r("6544"),h=r.n(p),f=r("7496"),m=r("40dc"),b=r("5bc1"),g=r("132d"),y=r("8860"),_=r("da13"),C=r("1baa"),I=r("34c3"),x=r("5d23"),S=r("f6c4"),V=r("f774"),k=r("2a7f"),w=Object(v["a"])(d,l,c,!1,null,null,null),M=w.exports;h()(w,{VApp:f["a"],VAppBar:m["a"],VAppBarNavIcon:b["a"],VIcon:g["a"],VList:y["a"],VListItem:_["a"],VListItemGroup:C["a"],VListItemIcon:I["a"],VListItemTitle:x["b"],VMain:S["a"],VNavigationDrawer:V["a"],VToolbarTitle:k["a"]});var E=r("8c4f"),R=function(){var e=this,t=e.$createElement,r=e._self._c||t;return r("v-container",[r("v-layout",[r("v-flex",{attrs:{xs12:""}},[r("h1",{staticClass:"display-1 mt-10 ml-10 mb-5"},[e._v("Welcome!")]),r("p",{staticClass:"ml-10"},[e._v(" The following tools are still to be tested/extended/etc. Use them with care. Feedback is greatly appreciated. ")]),r("v-row",{staticClass:"pl-10 pr-10"},[r("v-col",{attrs:{cols:"12",sm:"4",lg:"4"}},[r("v-hover",{scopedSlots:e._u([{key:"default",fn:function(t){var s=t.hover;return[r("v-card",{staticClass:"mx-auto transition-swing",attrs:{color:"grey lighten-5",elevation:s?4:2,to:{name:"NameChecker"}}},[r("v-card-text",{staticClass:"pt-6",staticStyle:{position:"relative"}},[r("h3",{staticClass:"display-1 font-weight-light blue--text mb-2"},[e._v(" Name Checker ")]),r("div",{staticClass:"font-weight-light title mb-2"},[e._v(" Takes a variant description as input and checks whether it is correct. ")])])],1)]}}])})],1),r("v-col",{attrs:{cols:"12",sm:"4",lg:"4"}},[r("v-hover",{scopedSlots:e._u([{key:"default",fn:function(t){var s=t.hover;return[r("v-card",{staticClass:"mx-auto transition-swing",attrs:{color:"grey lighten-4",elevation:s?4:2,to:{name:"PositionConverter"}}},[r("v-card-text",{staticClass:"pt-6",staticStyle:{position:"relative"}},[r("h3",{staticClass:"display-1 font-weight-light blue--text mb-2"},[e._v(" Position Converter ")]),r("div",{staticClass:"font-weight-light title mb-2"},[e._v(" Converts reference positions to selector orientated positions and vice versa. ")])])],1)]}}])})],1),r("v-col",{attrs:{cols:"12",sm:"4",lg:"4"}},[r("v-hover",{scopedSlots:e._u([{key:"default",fn:function(t){var s=t.hover;return[r("v-card",{staticClass:"mx-auto transition-swing",attrs:{color:"grey lighten-3",elevation:s?4:2,to:{name:"DescriptionExtractor"}}},[r("v-card-text",{staticClass:"pt-6",staticStyle:{position:"relative"}},[r("h3",{staticClass:"display-1 font-weight-light blue--text mb-2"},[e._v(" Description Extractor ")]),r("div",{staticClass:"font-weight-light title mb-2"},[e._v(" Generates the HGVS variant description from an observed and a reference sequence. ")])])],1)]}}])})],1)],1)],1)],1)],1)},O=[],q={name:"Home"},T=q,N=(r("e076"),r("b0af")),P=r("99d9"),$=r("62ad"),A=r("a523"),D=r("0e8f"),G=r("ce87"),z=r("a722"),L=r("0fd9"),j=Object(v["a"])(T,R,O,!1,null,"ed624c28",null),H=j.exports;h()(j,{VCard:N["a"],VCardText:P["a"],VCol:$["a"],VContainer:A["a"],VFlex:D["a"],VHover:G["a"],VLayout:z["a"],VRow:L["a"]});var F=function(){var e=this,t=e.$createElement,r=e._self._c||t;return r("v-container",[r("v-layout",[r("v-flex",{attrs:{xs12:""}},[r("h1",{staticClass:"display-1"},[e._v("About Page")]),r("p",[e._v(" Lorem ipsum dolor, sit amet consectetur adipisicing elit. Excepturi obcaecati tempora sunt debitis, minima deleniti ex inventore laboriosam at animi praesentium, quaerat corrupti molestiae recusandae corporis necessitatibus vitae, nam saepe? ")])])],1)],1)},B=[],U={},J=Object(v["a"])(U,F,B,!1,null,null,null),Y=J.exports;h()(J,{VContainer:A["a"],VFlex:D["a"],VLayout:z["a"]});var W,X=function(){var e=this,t=e.$createElement,r=e._self._c||t;return r("v-container",[r("v-layout",[r("v-flex",{attrs:{xs12:""}},[r("h1",{staticClass:"display-1 mt-10"},[e._v("Name Checker")]),r("p",[e._v(" The Name Checker takes the variant description as input and checks whether it is correct. ")]),r("v-sheet",{staticClass:"pa-10 mt-10",attrs:{elevation:"2"}},[r("v-text-field",{ref:"descriptionInput",attrs:{rules:e.rules,label:e.label,hint:e.hint,placeholder:e.placeholder,clearable:!0,autofocus:""},on:{keydown:function(t){return!t.type.indexOf("key")&&e._k(t.keyCode,"enter",13,t.key,"Enter")?null:e.nameCheck(t)}},model:{value:e.description,callback:function(t){e.description=t},expression:"description"}}),r("div",{staticClass:"examples-list"},[e._v(" Examples: "),e._l(e.examples,(function(t,s){return r("code",{key:s,on:{click:function(t){return t.preventDefault(),e.selectExample(s)}}},[r("span",{staticClass:"example-item"},[e._v(e._s(t))])])}))],2),r("v-btn",{ref:"nameCheck",staticClass:"mt-5",attrs:{disabled:!e.valid,color:"primary",to:{name:"NameChecker",params:{descriptionRouter:e.description}}}},[e._v(" Normalize ")]),r("v-overlay",{attrs:{absolute:!0,value:e.loadingOverlay}},[r("div",{staticClass:"text-center"},[r("v-progress-circular",{attrs:{size:50,indeterminate:""}})],1),r("div",{staticClass:"text-center"},[r("v-btn",{staticClass:"mt-5",on:{click:function(t){e.loadingOverlay=!1}}},[e._v(" Cancel ")])],1)])],1),e.summary?r("v-sheet",{staticClass:"pa-10 mt-10",attrs:{elevation:"2"}},[e.normalizedDescription?r("div",[r("h4",[e._v("Normalized Description")]),r("code",[r("span",{staticClass:"example-item"},[r("router-link",{attrs:{to:{name:"NameChecker",params:{descriptionRouter:e.normalizedDescription}}}},[e._v(e._s(e.normalizedDescription))])],1)])]):e._e()]):e._e()],1)],1)],1)},K=[],Q=(r("a4d3"),r("e01a"),i.a.CancelToken),Z={props:["descriptionRouter"],created:function(){this.descriptionRouter&&0!==this.descriptionRouter.length&&(this.description=this.descriptionRouter,this.nameCheck())},watch:{$route:function(){this.descriptionRouter&&0!==this.descriptionRouter.length&&(this.description=this.descriptionRouter,this.nameCheck())}},data:function(){return{description:null,descriptionModel:null,referenceModel:null,normalizedDescription:null,equivalentDescriptions:null,proteinDescriptions:null,visualize:null,sequence:null,summary:null,errors:null,loadingOverlay:!1,valid:!0,model:"",label:"HGVS Description",hint:"",placeholder:"",rules:[function(e){return!!e||"Required."}],examples:["NG_012337.1:g.7125G>T","LRG_24:g.5525_5532del"]}},methods:{selectExample:function(e){this.description=this.examples[e],this.$refs.descriptionInput.focus()},nameCheck:function(){var e=this;null!==this.description&&(this.loadingOverlay=!0,this.summary=null,this.errors=null,this.equivalentDescriptions=null,this.proteinDescriptions=null,this.visualize=null,this.descriptionModel=null,this.referenceModel=null,i.a.get("http://127.0.0.1:5000/api/name_check/"+this.description,{cancelToken:new Q((function(e){W=e}))}).then((function(t){t.data&&(e.summary=t.data,e.loadingOverlay=!1,e.normalizedDescription=t.data["normalized description"],t.data["errors"]&&(e.errors=t.data["errors"]),t.data["equivalent descriptions"]&&(e.equivalentDescriptions=t.data["equivalent descriptions"]),t.data["protein descriptions"]&&(e.proteinDescriptions=t.data["protein descriptions"]),t.data["visualize"]&&(e.visualize=t.data["visualize"]))})))},onShown:function(){this.$refs.cancel.focus()},onHidden:function(){},cancel:function(){W(),this.loadingOverlay=!1}}},ee=Z,te=(r("df01"),r("8336")),re=r("a797"),se=r("490a"),oe=r("8dd9"),ie=r("8654"),ne=Object(v["a"])(ee,X,K,!1,null,"88297008",null),ae=ne.exports;h()(ne,{VBtn:te["a"],VContainer:A["a"],VFlex:D["a"],VLayout:z["a"],VOverlay:re["a"],VProgressCircular:se["a"],VSheet:oe["a"],VTextField:ie["a"]});var le=function(){var e=this,t=e.$createElement,r=e._self._c||t;return r("v-container",[r("v-layout",[r("v-flex",{attrs:{xs12:""}},[r("h1",{staticClass:"display-1 mt-10"},[e._v("Position Converter")]),r("p",[e._v(" Converts reference positions to selector orientated positions and vice versa. ")]),r("v-sheet",{staticClass:"mt-10 pa-10",attrs:{elevation:"2"}},[r("v-form",{ref:"form",attrs:{"lazy-validation":e.lazy},model:{value:e.valid,callback:function(t){e.valid=t},expression:"valid"}},[r("v-row",[r("v-subheader",[e._v("From")])],1),r("v-divider"),r("v-row",[r("v-col",{attrs:{cols:"12",sm:"6",lg:"3"}},[r("v-text-field",{ref:"referenceId",attrs:{rules:e.rules,label:"Reference ID",hint:"E.g. NC_000001.11",clearable:!0,"error-messages":e.errorReferenceIdMessage,autofocus:""},model:{value:e.referenceId,callback:function(t){e.referenceId=t},expression:"referenceId"}})],1),r("v-col",{attrs:{cols:"12",sm:"6",lg:"3"}},[r("v-combobox",{ref:"fromSelectorId",attrs:{items:e.availableSelectors.selectors,label:"Selector ID",hint:"E.g. NM_001232.3","error-messages":e.errorSelectorIdMessage,clearable:!0},on:{click:function(t){return e.getAvailableSelectors()}},model:{value:e.fromSelectorId,callback:function(t){e.fromSelectorId=t},expression:"fromSelectorId"}})],1),r("v-col",{attrs:{cols:"12",sm:"6",lg:"3"}},[r("v-select",{attrs:{items:["Reference","Selector","g","c","n"],label:"Coordinate system"},model:{value:e.fromCoordinateSystem,callback:function(t){e.fromCoordinateSystem=t},expression:"fromCoordinateSystem"}})],1),r("v-col",{attrs:{cols:"12",sm:"6",lg:"3"}},[r("v-text-field",{ref:"position",attrs:{rules:e.rules,label:"Position",hint:"E.g. 100","error-messages":e.errorPositionMessage,clearable:!0},model:{value:e.position,callback:function(t){e.position=t},expression:"position"}})],1)],1),r("v-row",[r("v-subheader",[e._v("To")])],1),r("v-divider"),r("v-row",[r("v-col",[r("v-combobox",{attrs:{items:["Reference",this.availableSelectors.selectors],label:"Reference or selector ID"},model:{value:e.toSelectorId,callback:function(t){e.toSelectorId=t},expression:"toSelectorId"}})],1),r("v-col",{attrs:{cols:"12",sm:"6",lg:"3"}},[r("v-select",{attrs:{items:["Reference","Selector","g","c","n"],label:"Coordinate system"},model:{value:e.toCoordinateSystem,callback:function(t){e.toCoordinateSystem=t},expression:"toCoordinateSystem"}})],1)],1),r("v-row",{staticClass:"mt--10"},[r("v-col",[r("v-switch",{attrs:{label:"Include overlapping",color:"primary"},model:{value:e.includeOverlapping,callback:function(t){e.includeOverlapping=t},expression:"includeOverlapping"}})],1)],1)],1),r("v-row",[r("v-btn",{ref:"convert",staticClass:"mt-5",attrs:{color:"primary",disabled:!e.valid,to:{name:"PositionConverter",query:{referenceId:e.referenceId,fromSelectorId:e.fromSelectorId,fromCoordinateSystem:e.fromCoordinateSystem,position:e.position,toSelectorId:e.toSelectorId,toCoordinateSystem:e.toCoordinateSystem,includeOverlapping:e.includeOverlapping}}}},[e._v(" Convert ")]),r("v-spacer"),r("v-menu",{attrs:{transition:"slide-x-transition"},scopedSlots:e._u([{key:"activator",fn:function(t){var s=t.on,o=t.attrs;return[r("v-btn",e._g(e._b({staticClass:"mt-5",attrs:{color:"success"}},"v-btn",o,!1),s),[e._v(" Examples ")])]}}])},[r("v-list",e._l(e.examples,(function(t,s){return r("v-list-item",{key:s,attrs:{link:""}},[r("v-list-item-title",{attrs:{color:"success"},domProps:{textContent:e._s(t.item)},on:{click:function(r){return r.preventDefault(),e.setExample(t.fields)}}})],1)})),1)],1)],1)],1),e.errorMessages?r("v-alert",{staticClass:"mt-10",attrs:{border:"right",color:"red lighten-2",dark:"","colored-border":"",type:"error",elevation:"2",tile:""}},[e._v(" "+e._s(e.errorMessages)+" ")]):e._e(),e.summary?r("v-sheet",{staticClass:"pa-10 mt-10",attrs:{elevation:"2"}},[e._v(" "+e._s(e.summary)+" ")]):e._e(),e.summary?r("v-sheet",{staticClass:"pa-10 mt-10",attrs:{elevation:"2"}},[r("v-row",[r("v-col",[r("div",{staticClass:"overline mb-4"},[e._v("Reference Position")]),r("div",{},[r("span",{staticClass:"description"},[e._v(" "+e._s(e.summary.reference.id)+":"+e._s(e.summary.reference.coordinate_system)+"."+e._s(e.summary.reference.position)+" ")])])]),r("v-col",[r("div",{staticClass:"overline mb-4"},[e._v("Selector Position")]),r("div",[r("span",{staticClass:"description"},[e._v(" "+e._s(e.summary.reference.id)+"("+e._s(e.summary.selector.id)+"):"+e._s(e.summary.selector.coordinate_system)+"."+e._s(e.summary.selector.position)+" ")])])])],1),e.otherSelectors?r("v-row",[r("v-col",[r("div",{staticClass:"overline mb-4"},[e._v("Other Selectors")]),e.otherSelectors&&e.otherSelectors.length?e._e():r("div",[e._v(" No other selectors overlap with the provided position. ")]),e._l(e.otherSelectors,(function(t,s){return r("div",{key:s},[r("span",{staticClass:"description"},[e._v(" "+e._s(e.summary.reference.id)+"("+e._s(t.id)+"):"+e._s(t.coordinate_system)+"."+e._s(t.position)+" ")])])}))],2)],1):e._e()],1):e._e()],1)],1)],1)},ce=[],ue=r("b85c"),de={data:function(){return{valid:!0,lazy:!1,referenceId:"",fromSelectorId:"",fromCoordinateSystem:"",position:"",toSelectorId:"",toCoordinateSystem:"",includeOverlapping:!1,rules:[function(e){return!!e||"Required."}],summary:null,responseApi:null,errorMessages:null,errorReferenceId:null,errorReferenceIdMessage:null,errorSelectorId:null,errorSelectorIdMessage:null,errorPosition:null,errorPositionMessage:null,referenceErrors:null,otherSelectors:null,availableSelectors:{},examples:[{item:"LRG_24:g.100 -> t1",fields:{referenceId:"LRG_1",selectorId:"t1",position:"100",relativeTo:"Reference"}},{item:"NG_007485.1(NR_047542.1):n.274 -> NR_047542.1",fields:{referenceId:"NG_007485.1",selectorId:"NR_047542.1",position:"274",relativeTo:"Selector"}},{item:"NG_017013.2(NM_000546.5):c.274 -> NM_000546.5",fields:{referenceId:"NG_017013.2",selectorId:"NM_000546.5",position:"274",relativeTo:"Selector"}},{item:"NG_012337.1(NM_003002.2):c.274 -> NG_012337.1",fields:{referenceId:"NG_012337.1",selectorId:"NM_003002.2",position:"274",relativeTo:"Selector"}},{item:"NC_000001.11(NM_001232.3):c.274 -> NC_000001.11",fields:{referenceId:"NC_000001.11",selectorId:"NM_001232.3",position:"274",relativeTo:"Selector"}}]}},created:function(){this.run()},watch:{$route:function(){this.run()},referenceId:function(){this.availableSelectors.reference!==this.referenceId&&(this.availableSelectors={}),this.handleEretr(),this.handleEnoselector()},selectorId:function(){this.handleEnoselector()},position:function(){this.updatePositionErrorMessage()}},methods:{run:function(){this.$route.query.referenceId&&0!==this.$route.query.referenceId.length&&(this.referenceId=this.$route.query.referenceId),this.$route.query.selectorId&&0!==this.$route.query.selectorId.length&&(this.selectorId=this.$route.query.selectorId),this.$route.query.position&&0!==this.$route.query.position.length&&(this.position=this.$route.query.position),this.$route.query.relativeTo&&0!==this.$route.query.relativeTo.length&&(this.relativeTo=this.$route.query.relativeTo),this.$route.query.includeOverlapping&&(this.includeOverlapping=this.$route.query.includeOverlapping),this.$route.query.referenceId&&0!==this.$route.query.referenceId.length&&this.$route.query.selectorId&&0!==this.$route.query.selectorId.length&&this.$route.query.position&&0!==this.$route.query.position.length&&this.$route.query.relativeTo&&0!==this.$route.query.relativeTo.length&&this.positionConvert()},setExample:function(e){this.referenceId=e.referenceId,this.selectorId=e.selectorId,this.position=e.position,this.relativeTo=e.relativeTo},positionConvert:function(){var e=this;if(this.errorMessages=null,null!==this.referenceId&&null!==this.selectorId&&null!==this.position){this.summary=null,this.otherSelectors=null,this.responseApi=null;var t={reference_id:this.referenceId,selector_id:this.selectorId,position:this.position,relative_to:this.relativeTo,include_overlapping:this.includeOverlapping};i.a.get("http://127.0.0.1:5000/api/position_convert/",{params:t},{}).then((function(t){t.data&&(e.responseApi=t.data,e.responseHandler(t.data),e.loadingOverlay=!1)})).catch((function(t){t.response?e.errorMessages="Some response error.":t.request?e.errorMessages="Some request error.":e.errorMessages="Some error."}))}},responseHandler:function(e){"errors"in e?this.errorsHandler(e.errors):(this.summary=e,this.summary.other_selectors&&(this.otherSelectors=this.summary.other_selectors))},errorsHandler:function(e){var t,r=Object(ue["a"])(e);try{for(r.s();!(t=r.n()).done;){var s=t.value;"ERETR"===s.code?(this.errorMessages="Unable to retrieve reference "+this.referenceId,this.referenceErrors=["Unable to retrieve reference "+this.referenceId],this.errorReferenceId=this.referenceId,this.handleEretr()):"ENOSELECTOR"===s.code?(this.errorMessages="Selector "+this.selectorId+" not found in reference "+this.referenceId,this.errorSelectorId=this.selectorId,this.errorReferenceId=this.referenceId,this.handleEnoselector()):"ESYNTAX"===s.code?this.errorMessages="Position syntax error.":"ERANGELOCATION"===s.code?(this.errorMessages="",this.errorPosition=this.position,this.updatePositionErrorMessage("Range locations not supported.")):"EOUTOFBOUNDARY"===s.code?(this.errorPosition=this.position,this.updatePositionErrorMessage("Position out of sequence boundaries.")):s.code&&(this.errorMessages=s.code+" occurred.")}}catch(o){r.e(o)}finally{r.f()}},handleEretr:function(){this.errorReferenceId&&this.errorReferenceId===this.referenceId?this.errorReferenceIdMessage=this.referenceId+" not retrieved":this.errorReferenceIdMessage=null},handleEnoselector:function(){if(!this.errorReferenceId||this.errorReferenceId!==this.referenceId)return this.errorSelectorIdMessage=null,void(this.errorReferenceIdMessage=null);this.errorReferenceIdMessage=this.referenceId+" does not contain "+this.selectorId,this.errorSelectorId&&this.errorSelectorId===this.selectorId?this.errorSelectorIdMessage=this.selectorId+" not found in "+this.referenceId:(this.errorReferenceIdMessage=null,this.errorSelectorIdMessage=null)},updatePositionErrorMessage:function(e){this.errorPosition&&this.errorPosition===this.position?this.errorPositionMessage=e:this.errorPositionMessage=null},getAvailableSelectors:function(){var e=this;null!==this.referenceId&&0!==this.referenceId.length&&this.availableSelectors&&this.availableSelectors.reference!==this.referenceId&&(console.log("getAvailableSelectorsCalled"),i.a.get("http://127.0.0.1:5000/api/get_selectors/"+this.referenceId,{}).then((function(t){t.data&&(e.availableSelectors=t.data)})))}}},ve=de,pe=(r("dfb0"),r("0798")),he=r("2b5d"),fe=r("ce7e"),me=r("4bd4"),be=r("e449"),ge=r("b974"),ye=r("2fa4"),_e=r("e0c7"),Ce=r("b73d"),Ie=Object(v["a"])(ve,le,ce,!1,null,"0d85f774",null),xe=Ie.exports;h()(Ie,{VAlert:pe["a"],VBtn:te["a"],VCol:$["a"],VCombobox:he["a"],VContainer:A["a"],VDivider:fe["a"],VFlex:D["a"],VForm:me["a"],VLayout:z["a"],VList:y["a"],VListItem:_["a"],VListItemTitle:x["b"],VMenu:be["a"],VRow:L["a"],VSelect:ge["a"],VSheet:oe["a"],VSpacer:ye["a"],VSubheader:_e["a"],VSwitch:Ce["a"],VTextField:ie["a"]});var Se=function(){var e=this,t=e.$createElement,r=e._self._c||t;return r("v-container",[r("v-layout",[r("v-flex",{attrs:{xs12:""}},[r("h1",{staticClass:"display-1  mt-10"},[e._v("Description Extractor")]),r("p",[e._v(" Generates the HGVS variant description from a reference sequence and an observed sequence. ")]),r("v-sheet",{staticClass:"mt-10 pa-10",attrs:{elevation:"2"}},[r("v-form",{ref:"form",attrs:{"lazy-validation":e.lazy},model:{value:e.valid,callback:function(t){e.valid=t},expression:"valid"}},[r("v-row",[r("v-col",{attrs:{cols:"12"}},[r("v-text-field",{ref:"reference",attrs:{rules:e.rules,label:"Reference sequence",hint:"E.g. AATTTCCCGGG",clearable:!0,autofocus:""},model:{value:e.reference,callback:function(t){e.reference=t},expression:"reference"}})],1)],1),r("v-row",[r("v-col",{attrs:{cols:"12"}},[r("v-text-field",{ref:"observed",attrs:{rules:e.rules,label:"Observed Sequence",hint:"E.g. AATCCGG",clearable:!0},model:{value:e.observed,callback:function(t){e.observed=t},expression:"observed"}})],1)],1),r("v-row",[r("v-btn",{ref:"extract",staticClass:"mt-5",attrs:{color:"primary",disabled:!e.valid},on:{click:function(t){return t.preventDefault(),e.descriptionExtract()}}},[e._v(" Extract ")]),r("v-spacer"),r("v-btn",{staticClass:"mt-5",attrs:{color:"success"},on:{click:function(t){return t.preventDefault(),e.setExample()}}},[e._v(" Example ")])],1)],1)],1),e.summary?r("v-sheet",{staticClass:"pa-10 mt-10",attrs:{elevation:"2"}},[e._v(" "+e._s(e.summary)+" ")]):e._e()],1)],1)],1)},Ve=[],ke={data:function(){return{valid:!0,lazy:!1,reference:"",observed:"",position:"",rules:[function(e){return!!e||"Required."}],summary:null,responseApi:null,example:{reference:"AAA",observed:"ATA"}}},methods:{setExample:function(){this.reference=this.example.reference,this.observed=this.example.observed},descriptionExtract:function(){var e=this;if(null!==this.reference&&null!==this.observed){this.summary=null,this.responseApi=null;var t={reference:this.reference,observed:this.observed};i.a.get("http://127.0.0.1:5000/api/description_extract/",{params:t},{}).then((function(t){t.data&&(e.responseHandler(t.data),e.loadingOverlay=!1)})).catch((function(t){t.response?e.errorMessages="Some response error.":t.request?e.errorMessages="Some request error.":e.errorMessages="Some error."}))}},responseHandler:function(e){this.summary=e,"errors"in e?(console.log("some error"),this.errorsHandler(e.errors)):(console.log("we have some response"),this.summary=e)},errorsHandler:function(e){console.log("Some Error",e)}}},we=ke,Me=(r("cee8"),Object(v["a"])(we,Se,Ve,!1,null,"1def26ad",null)),Ee=Me.exports;h()(Me,{VBtn:te["a"],VCol:$["a"],VContainer:A["a"],VFlex:D["a"],VForm:me["a"],VLayout:z["a"],VRow:L["a"],VSheet:oe["a"],VSpacer:ye["a"],VTextField:ie["a"]}),s["a"].use(E["a"]);var Re=[{path:"/",name:"Home",component:H},{path:"/about",name:"About",component:Y},{path:"/namechecker/:descriptionRouter?",props:!0,name:"NameChecker",component:ae},{path:"/positionconverter",name:"PositionConverter",component:xe},{path:"/descriptionextractor",name:"DescriptionExtractor",component:Ee},{path:"*",name:"catchAll",component:H}],Oe=new E["a"]({mode:"history",base:"",routes:Re}),qe=Oe,Te=r("2f62");s["a"].use(Te["a"]);var Ne=new Te["a"].Store({state:{},mutations:{},actions:{},modules:{}}),Pe=r("f309");s["a"].use(Pe["a"]);var $e=new Pe["a"]({});s["a"].config.productionTip=!1,new s["a"]({router:qe,store:Ne,vuetify:$e,render:function(e){return e(M)}}).$mount("#app")},cee8:function(e,t,r){"use strict";var s=r("dcd3"),o=r.n(s);o.a},dcd3:function(e,t,r){},df01:function(e,t,r){"use strict";var s=r("e7d7"),o=r.n(s);o.a},dfb0:function(e,t,r){"use strict";var s=r("29bc"),o=r.n(s);o.a},e076:function(e,t,r){"use strict";var s=r("2f3f"),o=r.n(s);o.a},e7d7:function(e,t,r){}});
//# sourceMappingURL=app.b5a8918c.js.map