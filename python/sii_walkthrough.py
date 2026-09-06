"""教学图的计算：直接调用SII物理模块，保留同一基准场景及同一条随机记录。"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sii_unified as sii
from sii_observation import (source_coherence_spectrum, tracking_geometry,
                             simulate_array_photon_times)
from sii_performance import integer_align
from sii_validation import analytic_waveform_calibration


def theoretical_scene(root, instrument, observation, layout, source, figures, output):
    """先计算亮度、连续UV功率与飞秒HBT峰；还没有随机观测。"""
    root, figures, output = Path(root), Path(figures), Path(output)
    positions=layout[['east_m','north_m','up_m']].to_numpy()[:2]
    dec,lat=np.deg2rad([observation.source_dec_deg,observation.site_lat_deg])
    hour=.5  # 选观测窗口内有明显几何时差的时刻；全流程参数表明确记录。
    state=tracking_geometry(positions,hour,dec,lat)
    spectrum=source_coherence_spectrum(positions,hour,dec,lat,source,instrument,'single_disk')
    projected=sii.uvw_from_enu(positions[1]-positions[0],hour,dec,lat)
    baseline=np.hypot(projected[0],projected[1])
    theta=np.linspace(-.18,.18,361)
    x,y=np.meshgrid(theta,theta)
    sky=(x*x+y*y<=(source.primary_diameter_mas/2)**2).astype(float)
    sky/=sky.sum()
    lengths=np.linspace(0,1200,1201)
    visibility=sii.uniform_disk_visibility(lengths/instrument.wavelength_m,source.primary_diameter_mas)
    # 与主程序使用同一窄带及三条响应曲线；共同通光比例在归一场相关中抵消。
    wavelength=np.linspace(instrument.wavelength_nm-instrument.optical_width_nm/2,
                           instrument.wavelength_nm+instrument.optical_width_nm/2,4097)
    response=np.ones_like(wavelength)
    for filename in ['mirror_reflectivity_dm0113_13point_mean.csv','filter_transmission.csv','sipm_pde.csv']:
        xx,yy=sii._read_two_column_curve(root/'configs/efficiency'/filename)
        response*=np.interp(wavelength,xx,yy)
    density=response/wavelength
    lag_ps=np.linspace(-1.,1.,801)
    frequency=sii.C_M_S/(wavelength*1e-9)
    phase=np.exp(-2j*np.pi*(frequency-frequency.mean())[:,None]*lag_ps[None,:]*1e-12)
    mono_visibility=sii.uniform_disk_visibility(baseline/(wavelength*1e-9),source.primary_diameter_mas)
    field=np.trapezoid(density[:,None]*phase,wavelength,axis=0)/np.trapezoid(density,wavelength)
    field_pair=np.trapezoid((density*mono_visibility)[:,None]*phase,wavelength,axis=0)/np.trapezoid(density,wavelength)
    g2_zero=1+instrument.polarization_factor*abs(field)**2
    g2_pair=1+instrument.polarization_factor*abs(field_pair)**2
    fig,axes=plt.subplots(1,3,figsize=(13,3.7))
    axes[0].imshow(sky,origin='lower',extent=[theta[0],theta[-1],theta[0],theta[-1]],cmap='magma')
    axes[0].set(xlabel='East [mas]',ylabel='North [mas]',title=f'Sky: uniform disk {source.primary_diameter_mas:g} mas')
    axes[1].plot(lengths,visibility,label='V (real)');axes[1].plot(lengths,visibility**2,label='P = |V|^2')
    axes[1].axvline(baseline,ls=':',color='k',label='Tel.1 - Tel.2, H = 0.5 rad')
    axes[1].set(xlabel='Projected baseline [m]',ylabel='Visibility / squared visibility',title='Monochromatic theory, 400 nm')
    axes[1].legend(fontsize=7)
    axes[2].plot(lag_ps,g2_zero,label='Zero baseline');axes[2].plot(lag_ps,g2_pair,label='Tel.1 - Tel.2')
    axes[2].set(xlabel='Optical delay [ps]',ylabel='Ideal g^(2)',title='Before optical / electronic smearing')
    axes[2].legend(fontsize=8)
    fig.tight_layout();fig.savefig(figures/'09_theory_first.png');plt.show()
    pd.DataFrame({'lag_ps':lag_ps,'g2_zero':g2_zero,'g2_pair':g2_pair}).to_csv(output/'walkthrough_ideal_hbt.csv',index=False)
    values=dict(source_case='single_disk',diameter_mas=source.primary_diameter_mas,
        magnitude=source.ab_magnitude,wavelength_nm=instrument.wavelength_nm,width_nm=instrument.optical_width_nm,
        source_dec_deg=observation.source_dec_deg,site_lat_deg=observation.site_lat_deg,
        telescopes=len(layout),pair=[1,2],hour_angle_rad=hour,pair_projected_baseline_m=float(baseline),
        pair_power=float(spectrum['pair_visibility2'][0,1]),arrival_delay_ns=float(state['arrival_delays_ns'][1]),
        star_rate_hz=sii.detected_star_rate_hz(source.ab_magnitude,instrument),
        nsb_rate_hz=instrument.detected_nsb_rate_hz,coherence_area_s=instrument.coherence_area_s,
        sample_width_ns=instrument.sample_width_ns,raw_record_ns=24000.,raw_pair_scale=1.,
        hours=observation.hours_per_night,segment_s=observation.segment_s)
    (output/'walkthrough_parameters.json').write_text(json.dumps(values,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return dict(state=state,spectrum=spectrum,parameters=values)


def single_record(scene,instrument,figures,output,seed=20260910):
    """实际物理倍率1的同一对事件流供单镜图及后续两镜相关使用。"""
    params=scene['parameters'];rng=np.random.default_rng(seed)
    duration=params['raw_record_ns'];spectrum=scene['spectrum'];state=scene['state']
    template=sii.load_measured_spe_template(instrument.spe_template_path)
    events,metadata=simulate_array_photon_times(rng,duration,params['star_rate_hz'],params['nsb_rate_hz'],
        spectrum['coherence'],instrument,state['arrival_delays_ns'],pair_rate_scale=1.,
        padding_ns=float(np.max(abs(template[0]))),spectral_weights=spectrum['spectral_weights'])
    rendered=[sii.render_pe_waveform(rng,times,duration,instrument,template=template) for times in events]
    adc=np.stack([r['adc_mv'] for r in rendered]);times=rendered[0]['sample_time_ns']
    first_events=events[0][(events[0]>=0)&(events[0]<1000)]
    fig,axes=plt.subplots(3,1,figsize=(11,7),gridspec_kw={'height_ratios':[1,1.5,2]})
    axes[0].eventplot(first_events,lineoffsets=.5,linelengths=.7)
    axes[0].set(xlim=(0,1000),ylim=(0,1),yticks=[],xlabel='Time [ns]',title='Tel.1 detected event times: first 1 us')
    axes[1].plot(*template,color='C1');axes[1].set(xlabel='Time relative to one event [ns]',ylabel='SPE [mV]',title='Measured single-event response')
    for idx,label in [(0,'Tel.1'),(1,'Tel.2, original reception times')]:
        keep=times<1000;axes[2].plot(times[keep],adc[idx,keep],lw=.9,label=label)
    axes[2].set(xlabel='ADC time [ns]',ylabel='Sampled voltage [mV]',title='Same physical event streams: pulses overlap, sampled every 4 ns')
    axes[2].legend(fontsize=8);fig.tight_layout();fig.savefig(Path(figures)/'10_events_to_waveform.png');plt.show()
    np.savez_compressed(Path(output)/'walkthrough_record.npz',time_ns=times,adc_mv=adc,
        telescope_1_events_ns=events[0],telescope_2_events_ns=events[1],
        arrival_delays_ns=state['arrival_delays_ns'],pair_power=params['pair_power'],seed=seed)
    scene.update(adc=adc,time_ns=times,events=events,seed=seed)
    return pd.DataFrame([dict(telescope=i+1,events_inside_record=int(np.sum((ev>=0)&(ev<duration))),
        voltage_mean_mv=float(adc[i].mean()),voltage_rms_mv=float(adc[i].std())) for i,ev in enumerate(events)])


def pair_correlation(scene,instrument,figures,output):
    """同一记录的相关函数，不因物理峰被噪声淹没而放大注入。"""
    adc=scene['adc'];dt=instrument.sample_width_ns;delay=scene['state']['arrival_delays_ns']
    aligned,residual,exposure=integer_align(adc,delay,dt)
    lag_raw,raw=sii.waveform_cross_correlation(adc[0],adc[1],dt,400.)
    lag,observed=sii.waveform_cross_correlation(*aligned,dt,200.)
    calibration,_=analytic_waveform_calibration(instrument,block_duration_ns=exposure*1e9,residual_delay_ns=float(residual[0]))
    power=scene['parameters']['pair_power']
    expected=power*calibration.peak_per_visibility2
    estimate,sigma=sii.estimate_visibility2_gls(observed[None,:],calibration)
    fig,axes=plt.subplots(1,3,figsize=(13,3.6))
    axes[0].plot(lag_raw,raw,lw=1)
    axes[0].axvline(-delay[1],ls=':',color='C1',label='Expected geometric peak position')
    axes[0].set(xlabel='Lag [ns]',ylabel='Voltage correlation C',title='Before integer delay alignment');axes[0].legend(fontsize=7)
    axes[1].plot(lag,observed,label='One physical 24 us record');axes[1].plot(lag,expected,label='Expected signal')
    axes[1].set(xlabel='Lag after alignment [ns]',ylabel='Voltage correlation C',title='Physical signal is below short-record noise');axes[1].legend(fontsize=7)
    for duration in [1200.,21600.]:
        predicted_sigma=sii.waveform_gls_weights(calibration,duration)[1]
        axes[2].errorbar(duration/3600,power,yerr=predicted_sigma,fmt='o',label=f'{duration/3600:g} h: expected +/- 1 sigma')
    axes[2].set(xlabel='Integration [h]',ylabel='Estimated squared visibility',title='Fixed-baseline precision forecast');axes[2].legend(fontsize=7)
    fig.tight_layout();fig.savefig(Path(figures)/'11_pair_observation.png');plt.show()
    fig,axis=plt.subplots(figsize=(7,3))
    axis.plot(lag,expected*1e6);axis.set(xlabel='Lag after alignment [ns]',ylabel='Expected C [1e-6]',title='Expected broadened peak only; separate vertical scale')
    fig.tight_layout();fig.savefig(Path(figures)/'12_expected_correlation.png');plt.show()
    table=pd.DataFrame({'lag_ns':lag,'measured_C':observed,'expected_C':expected})
    table.to_csv(Path(output)/'walkthrough_pair_correlation.csv',index=False)
    result=dict(true_power=power,estimated_power=float(estimate[0]),sigma_short=float(sigma),
        effective_record_s=exposure,residual_delay_ns=float(residual[0]),
        expected_peak=float(expected.max()),sigma_1200s=sii.waveform_gls_weights(calibration,1200.)[1])
    (Path(output)/'walkthrough_pair_result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


def array_and_model_views(pipelines,instrument,observation,figures,output):
    """连续理论、实际采样和随机测量分开画；各模型明确列出参数变化。"""
    frame=pipelines['single_disk'].measurements
    limit=max(abs(frame.u_lambda).max(),abs(frame.v_lambda).max())/1e6*1.03
    axis=np.linspace(-limit,limit,321);uu,vv=np.meshgrid(axis*1e6,axis*1e6)
    theoretical=sii.uniform_disk_visibility(np.hypot(uu,vv),.16)**2
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    image=axes[0].imshow(theoretical,origin='lower',extent=[-limit,limit,-limit,limit],vmin=0,vmax=1,cmap='viridis')
    axes[0].set(title='Continuous theory: 0.16 mas disk')
    for ax,column,title in zip(axes[1:],['visibility2_true','visibility2_measured'],['Sampled, time + band average','Simulated measurements']):
        image=ax.scatter(frame.u_lambda/1e6,frame.v_lambda/1e6,c=frame[column],s=2,vmin=0,vmax=1,cmap='viridis')
        ax.set(title=title,xlim=(-limit,limit),ylim=(-limit,limit))
    for ax in axes:
        ax.set(xlabel='u [million wavelengths]',ylabel='v [million wavelengths]');ax.set_aspect('equal')
    fig.colorbar(image,ax=axes,label='Squared visibility P (display range 0 to 1)',fraction=.025)
    fig.savefig(Path(figures)/'13_disk_uv_observation.png',bbox_inches='tight');plt.show()
    frame.drop(columns=['uv_samples_u','uv_samples_v','uv_samples_weight']).to_csv(Path(output)/'walkthrough_disk_measurements.csv',index=False)


def model_views(pipelines,instrument,figures):
    """固定仪器和阵列，只改变明确列出的源参数。"""
    fig,axes=plt.subplots(1,4,figsize=(14,3.6))
    for ax,(case,pipeline) in zip(axes,pipelines.items()):
        data=pipeline.measurements
        ax.scatter(data.u_lambda/1e6,data.v_lambda/1e6,c=data.visibility2_true,s=2,vmin=0,vmax=1,cmap='viridis')
        ax.set(title=case,xlabel='u [million wavelengths]',ylabel='v [million wavelengths]');ax.set_aspect('equal')
    fig.tight_layout();fig.savefig(Path(figures)/'14_models_on_same_array.png');plt.show()
    lengths=np.linspace(0,1200,1201)
    fig,axes=plt.subplots(1,2,figsize=(10,3.6))
    for diameter in [.08,.16,.32]:
        axes[0].plot(lengths,sii.uniform_disk_visibility(lengths/instrument.wavelength_m,diameter)**2,label=f'{diameter:g} mas')
    for case,pipeline in pipelines.items():
        data=pipeline.measurements
        axes[1].hist(data.visibility2_true,bins=np.linspace(0,1,31),histtype='step',label=case)
    axes[0].set(xlabel='Projected baseline [m]',ylabel='Theoretical P',title='Only change the disk diameter')
    axes[1].set(xlabel='Time + band averaged P',ylabel='Number of sampled measurements',title='Same geometry, different source structures')
    for ax in axes: ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(Path(figures)/'15_parameter_effects.png');plt.show()
